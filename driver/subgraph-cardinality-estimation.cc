/*
#include <iostream>
#include <set>
#include "Base/Metrics.h"
#include "Base/Timer.h"
#include "DataStructure/Graph.h"
#include "SpecialSubgraphs/SmallCycle.h"
#include "SubgraphMatching/DataGraph.h"
#include "SubgraphMatching/PatternGraph.h"
#include "SubgraphMatching/CandidateSpace.h"
#include "SubgraphMatching/CandidateFilter.h"
#include "SubgraphCounting/Option.h"
#include "SubgraphCounting/CardinalityEstimation.h"

using namespace std;
using namespace GraphLib;

std::set<std::string> scientific_type_results = {"#CandTree"};
std::set<std::string> double_type_results = {
    "Truth", "Est", "logQError", "CSBuildTime", "TreeCountTime", "TreeSampleTime", "GraphSampleTime", "QueryTime"
};
std::set<std::string> longlong_type_results = {};
std::vector<std::string> print_order = {
    "#CSVertex", "#CSEdge", "#CandTree", "#TreeTrials", "#TreeSuccess","Truth", "Est", "logQError",
    "CSBuildTime", "TreeCountTime", "TreeSampleTime", "GraphSampleTime", "QueryTime"
};
std::vector<dict> results;
std::string query_path;
Timer timer;
std::vector<PatternGraph*> pattern_graphs;
std::deque<std::string> query_names;
std::unordered_map<std::string, double> true_cnt;
double total_time = 0.0;

void read_ans(const std::string& dataset) {
    std::string ans_file_name = query_path;
    cout << ans_file_name << endl;
    std::ifstream ans_in(ans_file_name);
    while (!ans_in.eof()) {
        std::string name, t, c;
        ans_in >> name >> t >> c;
        if (name.empty() || c.empty()) continue;
        name = "../dataset/"+dataset+"/query_graph/"+name;
        true_cnt[name] = stod(c);
        query_names.push_back(name);
    }
}


void read_filter_option(const std::string& opt, const std::string &filter, CardinalityEstimation::CardEstOption& option) {
    if (opt.substr(2) == "STRUCTURE") {
        if (filter == "X")
            option.structure_filter = SubgraphMatching::NO_STRUCTURE_FILTER;
        else if (filter == "3")
            option.structure_filter = SubgraphMatching::TRIANGLE_SAFETY;
        else if (filter == "4")
            option.structure_filter = SubgraphMatching::FOURCYCLE_SAFETY;
    }
}

int32_t main(int argc, char *argv[]) {
    std::string dataset = "wordnet";
    CardinalityEstimation::CardEstOption opt;

    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'd':
                    dataset = argv[i + 1];
                    break;
                case 'q':
                    query_path = argv[i + 1];
                    break;
                case 'K':
                    opt.ub_initial = atoi(argv[i + 1]);
                    break;
                case '-':
                    read_filter_option(std::string(argv[i]), std::string(argv[i+1]), opt);
                    break;
                default:
                    break;
            }
        }
    }

    if (query_path.empty()) {
        query_path = "../dataset/"+dataset+"/"+dataset+"_ans.txt";
    }

    std::string data_path = "../dataset/"+dataset+"/"+dataset+".graph";
    read_ans(dataset);
    DataGraph D;
    D.LoadLabeledGraph(data_path);
    D.Preprocess();
    opt.MAX_QUERY_VERTEX = 12;
    opt.MAX_QUERY_EDGE = 4;
    pattern_graphs.resize(query_names.size());
    for (int i = 0; i < query_names.size(); i++) {
        std::string query_name = query_names[i];

        //动态分配内存
        pattern_graphs[i] = new PatternGraph();

        pattern_graphs[i]->LoadLabeledGraph(query_name);
        pattern_graphs[i]->ProcessPattern(D);
        pattern_graphs[i]->EnumerateLocalTriangles();
        pattern_graphs[i]->EnumerateLocalFourCycles();
        opt.MAX_QUERY_VERTEX = std::max(opt.MAX_QUERY_VERTEX, pattern_graphs[i]->GetNumVertices());
        opt.MAX_QUERY_EDGE = std::max(opt.MAX_QUERY_EDGE, pattern_graphs[i]->GetNumEdges());
    }

    if (opt.structure_filter >= SubgraphMatching::FOURCYCLE_SAFETY) {
        D.EnumerateLocalFourCycles();
    }
    if (opt.structure_filter >= SubgraphMatching::TRIANGLE_SAFETY) {
        D.EnumerateLocalTriangles();
    }
    CardinalityEstimation::FaSTestCardinalityEstimation estimator(&D, opt);

    for (int i = 0; i < pattern_graphs.size(); i++) {
        PatternGraph* P = pattern_graphs[i];
        std::string query_name = query_names[i];
        std::cout << "Start Processing " << query_name << std::endl;
        double est = estimator.EstimateEmbeddings(P);
        dict query_result = estimator.GetResult();
        query_result["Est"] = est;
        if (true_cnt.find(query_name)!= true_cnt.end()) {
            query_result["Truth"] = std::any(true_cnt[query_name]*1.0);
            query_result["logQError"] = std::any(logQError(true_cnt[query_name]*1.0, est));
        }
        for (auto &key : print_order) {
            if (query_result.find(key) == query_result.end()) continue;
            std::any value = query_result[key];
            if (double_type_results.find(key) != double_type_results.end())
                fprintf(stdout, "  [Result] %-20s: %.04lf\n", key.c_str(), std::any_cast<double>(value));
            else if (scientific_type_results.find(key)!= scientific_type_results.end())
                fprintf(stdout, "  [Result] %-20s: %.04g\n", key.c_str(), std::any_cast<double>(value));
            else if (longlong_type_results.find(key)!= longlong_type_results.end())
                fprintf(stdout, "  [Result] %-20s: %lld\n", key.c_str(), std::any_cast<long long>(value));
            else
                fprintf(stdout, "  [Result] %-20s: %d\n", key.c_str(), std::any_cast<int>(value));
        }
        cout << query_name << " Finished!\n";
        fflush(stdout);
        results.push_back(query_result);
    }

    //增加：释放动态分配的 PatternGraph 内存
    for (int i = 0; i < pattern_graphs.size(); i++) {
        delete pattern_graphs[i];  // 释放内存
    }

    std::function<double(double)> absolute_value = [](double x) {return std::abs(x);};
    cout << std::fixed << std::setprecision(2) << "Total Time: " << Total(results, "QueryTime") << "ms\n";
    cout << std::fixed << std::setprecision(2) <<
        "Average Abs log Q-Error: " << Average(results, "logQError", absolute_value) << endl;
}
*/

//删除比对部分代码
/********************** 头文件 ***********************/
#include <iostream>
#include <set>
#include <future>
#include <chrono>
#include <optional>
#include "Base/Metrics.h"
#include "Base/Timer.h"
#include "DataStructure/Graph.h"
#include "SpecialSubgraphs/SmallCycle.h"
#include "SubgraphMatching/DataGraph.h"
#include "SubgraphMatching/PatternGraph.h"
#include "SubgraphMatching/CandidateSpace.h"
#include "SubgraphMatching/CandidateFilter.h"
#include "SubgraphCounting/Option.h"
#include "SubgraphCounting/CardinalityEstimation.h"

/************************** 命名空间 *************************/
using namespace std;
using namespace GraphLib;

/**************************** 全局变量 *************************/
std::set<std::string> scientific_type_results = {"#CandTree"};
std::set<std::string> double_type_results = {
     "Est", "CSBuildTime", "TreeCountTime", "TreeSampleTime", "GraphSampleTime", "QueryTime"
};
std::set<std::string> longlong_type_results = {};
std::vector<std::string> print_order = {
    "#CSVertex", "#CSEdge", "#CandTree", "#TreeTrials", "#TreeSuccess", "Est",
    "CSBuildTime", "TreeCountTime", "TreeSampleTime", "GraphSampleTime", "QueryTime"
};

std::vector<dict> results;
std::string query_path;
Timer timer;
std::vector<PatternGraph*> pattern_graphs;
std::deque<std::string> query_names;
double total_time = 0.0;

struct ProcessResult {
    dict query_result;
    double est;
};

void read_ans(const std::string& dataset) {
    std::string ans_file_name = query_path;
    cout << "Reading from: " << ans_file_name << endl;
    ifstream ans_in(ans_file_name);
    if (!ans_in) {
        cout << "Error: Unable to open file " << ans_file_name << endl;
        return;
    }
    while (!ans_in.eof()) {
        string name;
        ans_in >> name;
        if (name.empty()) continue;
        name = "../dataset/" + dataset + "/query_graph/" + name;
        query_names.push_back(name);
    }
}

void read_filter_option(const string& opt, const string &filter, CardinalityEstimation::CardEstOption& option) {
    if (opt.substr(2) == "STRUCTURE") {
        if (filter == "X") option.structure_filter = SubgraphMatching::NO_STRUCTURE_FILTER;
        else if (filter == "3") option.structure_filter = SubgraphMatching::TRIANGLE_SAFETY;
        else if (filter == "4") option.structure_filter = SubgraphMatching::FOURCYCLE_SAFETY;
    }
}

int32_t main(int argc, char *argv[]) {
    string dataset = "wordnet";
    CardinalityEstimation::CardEstOption opt;

    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'd': dataset = argv[i + 1]; break;
                case 'q': query_path = argv[i + 1]; break;
                case 'K': opt.ub_initial = atoi(argv[i + 1]); break;
                case '-': read_filter_option(string(argv[i]), string(argv[i+1]), opt); break;
            }
        }
    }

    if (query_path.empty()) {
        query_path = "../dataset/"+dataset+"/"+dataset+"_ans.txt";
    }

    string data_path = "../dataset/"+dataset+"/"+dataset+".graph";
    read_ans(dataset);
    DataGraph D;
    D.LoadLabeledGraph(data_path);
    D.Preprocess();
    opt.MAX_QUERY_VERTEX = 12;
    opt.MAX_QUERY_EDGE = 4;

    pattern_graphs.resize(query_names.size());
    for (size_t i = 0; i < query_names.size(); i++) {
        string query_name = query_names[i];
        pattern_graphs[i] = new PatternGraph();
        pattern_graphs[i]->LoadLabeledGraph(query_name);
        pattern_graphs[i]->ProcessPattern(D);
        pattern_graphs[i]->EnumerateLocalTriangles();
        pattern_graphs[i]->EnumerateLocalFourCycles();
        opt.MAX_QUERY_VERTEX = max(opt.MAX_QUERY_VERTEX, pattern_graphs[i]->GetNumVertices());
        opt.MAX_QUERY_EDGE = max(opt.MAX_QUERY_EDGE, pattern_graphs[i]->GetNumEdges());
    }

    if (opt.structure_filter >= SubgraphMatching::FOURCYCLE_SAFETY) D.EnumerateLocalFourCycles();
    if (opt.structure_filter >= SubgraphMatching::TRIANGLE_SAFETY) D.EnumerateLocalTriangles();

    /*CardinalityEstimation::FaSTestCardinalityEstimation estimator(&D, opt);

    for (size_t i = 0; i < pattern_graphs.size(); i++) {
        PatternGraph* P = pattern_graphs[i];
        string query_name = query_names[i];
        cout << "Start Processing " << query_name << endl;

        auto future = async(launch::async, [&D, opt, P]() -> optional<ProcessResult> {
            try {
                CardinalityEstimation::FaSTestCardinalityEstimation local_estimator(&D, opt);
                double est = local_estimator.EstimateEmbeddings(P);
                dict query_result = local_estimator.GetResult();
                query_result["Est"] = est;
                return ProcessResult{query_result, est};
            } catch (...) {
                return nullopt;
            }
        });

        future_status status = future.wait_for(chrono::minutes(5));

        if (status == future_status::timeout) {
            cerr << "Error: Processing " << query_name 
                 << " timed out after 5 minutes. Skipping to next file." << endl;
            continue;
        }

        auto result_opt = future.get();
        if (!result_opt.has_value()) {
            cerr << "Error processing " << query_name << ". Skipping." << endl;
            continue;
        }

        auto& [query_result, est] = result_opt.value();
        for (auto &key : print_order) {
            if (query_result.find(key) == query_result.end()) continue;
            any value = query_result[key];
            if (double_type_results.count(key)) {
                printf("  [Result] %-20s: %.04lf\n", key.c_str(), any_cast<double>(value));
            } else if (scientific_type_results.count(key)) {
                printf("  [Result] %-20s: %.04g\n", key.c_str(), any_cast<double>(value));
            } else if (longlong_type_results.count(key)) {
                printf("  [Result] %-20s: %lld\n", key.c_str(), any_cast<long long>(value));
            } else {
                printf("  [Result] %-20s: %d\n", key.c_str(), any_cast<int>(value));
            }
        }
        cout << query_name << " Finished!\n";
        results.push_back(query_result);
    }

    cout << fixed << setprecision(2) << "Total Time: " << Total(results, "QueryTime") << "ms\n";
*/

/*
* @param
* D : 
* opt : 
*/
CardinalityEstimation::FaSTestCardinalityEstimation estimator(&D, opt);

for (size_t i = 0; i < pattern_graphs.size(); i++) {
    PatternGraph* P = pattern_graphs[i];
    string query_name = query_names[i];
    cout << "Start Processing " << query_name << endl;

    // ================== 修改点：移除异步和多线程代码 ==================
    try {
        // 直接创建 estimator 并执行估算（原异步任务中的代码）
        CardinalityEstimation::FaSTestCardinalityEstimation local_estimator(&D, opt);
        double est = local_estimator.EstimateEmbeddings(P);
        dict query_result = local_estimator.GetResult();
        query_result["Est"] = est;

        // 输出结果（与原逻辑一致）
        for (auto &key : print_order) {
            if (query_result.find(key) == query_result.end()) continue;
            any value = query_result[key];
            if (double_type_results.count(key)) {
                printf("  [Result] %-20s: %.04lf\n", key.c_str(), any_cast<double>(value));
            } else if (scientific_type_results.count(key)) {
                printf("  [Result] %-20s: %.04g\n", key.c_str(), any_cast<double>(value));
            } else if (longlong_type_results.count(key)) {
                printf("  [Result] %-20s: %lld\n", key.c_str(), any_cast<long long>(value));
            } else {
                printf("  [Result] %-20s: %d\n", key.c_str(), any_cast<int>(value));
            }
        }
        cout << query_name << " Finished!\n";
        results.push_back(query_result);
    } catch (...) {
        cerr << "Error processing " << query_name << ". Skipping." << endl;
        continue;
    }
    // ==============================================================
}

cout << fixed << setprecision(2) << "Total Time: " << Total(results, "QueryTime") << "ms\n";

// 添加释放 PatternGraph 内存的代码（修复内存泄漏）
for (auto* pg : pattern_graphs) {
    delete pg;
}
pattern_graphs.clear();

return 0;
}
