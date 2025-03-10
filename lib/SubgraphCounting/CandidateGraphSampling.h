#pragma once

namespace GraphLib {
    namespace CardinalityEstimation {
        class CandidateGraphSampler {
            CardEstOption opt;
            CandidateSpace *CS;
            DataGraph *data_;
            PatternGraph *query_;
            dict info;
            int root;
            int min_cand;
            bool* seen;
            int **local_candidates, *local_candidate_size;

            int num_embeddings = 0;

        public:
            dict GetInfo() {return info;}
            CandidateGraphSampler(DataGraph *data, CardEstOption opt_){
                data_ = data;
                seen = new bool[data->GetNumVertices()];
                opt = opt_;
                local_candidates = new int*[opt.MAX_QUERY_VERTEX];
                for (int i = 0; i < opt.MAX_QUERY_VERTEX; i++) {
                    local_candidates[i] = new int[data->GetNumVertices()];
                }
                local_candidate_size = new int[opt.MAX_QUERY_VERTEX];
            };
            ~CandidateGraphSampler(){
                for (int i = 0; i < opt.MAX_QUERY_VERTEX; i++) {
                    delete[] local_candidates[i];
                }
                delete[] local_candidates;
                delete[] local_candidate_size;
                delete[] seen;
            };

            void Preprocess(PatternGraph *query, CandidateSpace *cs) {
                info.clear();
                query_ = query;
                CS = cs;
                min_cand = (CS->GetNumCSVertex() > 1e4 || query->GetNumVertices() > 20) ? 2 : 4;
                sample.resize(query_->GetNumVertices(), -1);
                memset(seen, 0, data_->GetNumVertices());
                memset(local_candidate_size, 0, query->GetNumVertices());
            }
            std::vector<int> sample;
            std::vector<std::pair<std::vector<int>::iterator, std::vector<int>::iterator>> iterators;
            std::vector<int> root_candidates_;

            double Estimate(int ub_initial);

            std::pair<double, int> StratifiedSampling(int vertex_id, int ub, double w);

            int ChooseExtendableVertex();

            void BuildExtendableCandidates(int u);

            void Intersection(int index);
        };

        int printcnt = 0;
        double CandidateGraphSampler::Estimate(int ub_initial) {
            printcnt = 0;
            Timer timer; timer.Start();
            std::vector <int> num_cands(query_->GetNumVertices());
            for (int i = 0; i < query_->GetNumVertices(); i++) {
                num_cands[i] = CS->GetCandidateSetSize(i);
            }
            root = std::min_element(num_cands.begin(), num_cands.end()) - num_cands.begin();
            root_candidates_ = CS->GetCandidates(root);
            std::shuffle(root_candidates_.begin(), root_candidates_.end(), gen);
            double est = 0.0;
            int ub = ub_initial;
            int num_root_samples = root_candidates_.size();
            int used_samples = 0;
            for (int i = 0; i < num_root_samples; i++) {
                std::fill(sample.begin(), sample.end(), -1);
                memset(local_candidate_size, 0, query_->GetNumVertices());
                sample[root] = (i % root_candidates_.size());
                int num_sample_use = (ub - used_samples) / (std::max(num_root_samples - i, 1));
                auto sampling_result = StratifiedSampling(1, num_sample_use, 1.0 * root_candidates_.size());
                est += sampling_result.first;
                used_samples += sampling_result.second;
            }
            est /= num_root_samples;
            timer.Stop();
            info["GraphSampleTime"] = timer.GetTime();
            return est;
        }

        std::pair<double, int> CandidateGraphSampler::StratifiedSampling(int vertex_id, int ub, double w) {
            int u = ChooseExtendableVertex();
            BuildExtendableCandidates(u);

            if (local_candidate_size[u] == 0) {
                return {0, 1};
            }
            for (int i = 0; i < query_->GetNumVertices(); i++) {
                if (sample[i] == -1) continue;
                seen[CS->GetCandidate(i, sample[i])] = true;
            }

            for (int i = 0; i < local_candidate_size[u]; ++i) {
                if (seen[CS->GetCandidate(u, local_candidates[u][i])]) {
                    local_candidates[u][i] = local_candidates[u][local_candidate_size[u]-1];
                    local_candidate_size[u]--;
                    i--;
                }
            }

            for (int i = 0; i < query_->GetNumVertices(); i++) {
                if (sample[i] == -1) continue;
                seen[CS->GetCandidate(i, sample[i])] = false;
            }

            if (local_candidate_size[u] == 0) {
                local_candidate_size[u] = 0;
                sample[u] = -1;
                return {0, 1};
            }
            if (vertex_id == query_->GetNumVertices()-1) {
                double return_value = local_candidate_size[u] * 1.0;
                local_candidate_size[u] = 0;
                sample[u] = -1;
                return {w * return_value, 1};
            }

            int sample_space_size = local_candidate_size[u];

            int num_strata = std::min(std::max((int)ceil(sample_space_size * opt.strata_ratio), min_cand), ub);
            num_strata = std::min(num_strata, sample_space_size);

            int num_used = 0;
            double est = 0.0;
            if (num_strata == -1) {
                sample[u] = local_candidates[u][gen()%local_candidate_size[u]];
                std::tie(est, num_used) = StratifiedSampling(vertex_id + 1, ub, w * sample_space_size);
                sample[u] = -1;
                local_candidate_size[u] = 0;
                return {est, num_used};
            }
            else {
                int i = 0;
                while (num_used < ub and local_candidate_size[u] > 0) {
                    int idx = gen()%local_candidate_size[u];
                    sample[u] = local_candidates[u][idx];
                    int num_next_samples = (ub - num_used) / std::max(num_strata - i, 1);
                    if (num_next_samples == 0) num_next_samples = (ub - num_used);
                    double est_; int num_used_;
                    std::tie(est_, num_used_) = StratifiedSampling(vertex_id + 1, num_next_samples, w * sample_space_size);
                    est += est_;
                    num_used += num_used_;
                    local_candidates[u][idx] = local_candidates[u][local_candidate_size[u]-1];
                    local_candidate_size[u]--;
                    sample[u] = -1;
                    i++;
                    if (i == num_strata) break;
                }
                sample[u] = -1;
                local_candidate_size[u] = 0;
                return {est / i, num_used};
            }
        }

        int CandidateGraphSampler::ChooseExtendableVertex() {
            int u = -1;
            int max_open_neighbors = 0;
            int min_nbr_cnt = 1e9;
            for (int i = 0; i < query_->GetNumVertices(); i++) {
                if (sample[i] != -1) continue;
                int nbr_cnt = 1e9;
                int open_neighbors = 0;
                for (int q_nbr : query_->GetNeighbors(i)) {
                    if (sample[q_nbr] != -1) {
                        open_neighbors++;
                        int num_nbr = CS->GetCandidateNeighbors(q_nbr, sample[q_nbr], i).size();
                        if (num_nbr < nbr_cnt) {
                            nbr_cnt = num_nbr;
                        }
                    }
                }
                if (open_neighbors > max_open_neighbors) {
                    max_open_neighbors = open_neighbors;
                    min_nbr_cnt = nbr_cnt;
                    u = i;
                }
                else if (open_neighbors == max_open_neighbors) {
                    if (nbr_cnt < min_nbr_cnt) {
                        min_nbr_cnt = nbr_cnt;
                        u = i;
                    }
                }
            }
            return u;
        }

        void CandidateGraphSampler::BuildExtendableCandidates(int u) {
            local_candidate_size[u] = 0;
            iterators.clear();
            for (int q_nbr : query_->GetNeighbors(u)) {
                if (sample[q_nbr] == -1) continue;
                auto &candidate_neighbors = CS->GetCandidateNeighbors(q_nbr, sample[q_nbr], u);
                iterators.emplace_back(candidate_neighbors.begin(), candidate_neighbors.end());
            }
            std::sort(iterators.begin(), iterators.end(), [](auto &a, auto &b) -> bool {
                return a.second - a.first < b.second - b.first;
            });
            Intersection(u);
        }

        void CandidateGraphSampler::Intersection(int index) {
            if (local_candidate_size[index] > 0) return;
            int num_vectors = iterators.size();
            if (num_vectors == 1) {
                while (iterators[0].first != iterators[0].second) {
                    local_candidates[index][local_candidate_size[index]++] = (*iterators[0].first);
                    ++iterators[0].first;
                }
                return;
            }
            while (iterators[0].first != iterators[0].second) {
                int target = *iterators[0].first;
                for (int i = 1; i < num_vectors; i++) {
                    while (iterators[i].first != iterators[i].second) {
                        if (*iterators[i].first < target) {
                            ++iterators[i].first;
                        }
                        else if (*iterators[i].first > target) {
                            goto nxt_target;
                        }
                        else break;
                    }
                    if (iterators[i].first == iterators[i].second) return;
                }
                local_candidates[index][local_candidate_size[index]++] = target;
                nxt_target:
                ++iterators[0].first;
            }
        }

    }
}