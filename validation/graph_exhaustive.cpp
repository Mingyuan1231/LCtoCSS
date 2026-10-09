#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <vector>

// Independent comparison: affine recognition versus full labeled LC orbits.
// Orbit enumeration is independent of the affine solver.
using Mask = std::uint32_t;
using Row = std::uint64_t;
using Adj = std::array<Mask, 7>;

bool solve(std::vector<Row> equations, int d, Row &solution) {
    std::array<Row, 15> pivot{};
    const Row coefficients = (Row(1) << d) - 1;
    for (Row line : equations) {
        for (int c = 0; c < d; ++c) {
            if (!(line & (Row(1) << c))) continue;
            if (pivot[c]) line ^= pivot[c];
            else { pivot[c] = line; line = 0; break; }
        }
        if (!(line & coefficients) && (line & (Row(1) << d))) return false;
    }
    solution = 0;
    for (int c = d - 1; c >= 0; --c) {
        if (!pivot[c]) continue;
        const int value = int((pivot[c] >> d) & 1) ^
                          __builtin_parityll(pivot[c] & solution);
        if (value) solution |= Row(1) << c;
    }
    return true;
}

Adj adjacency(Mask graph, const std::vector<std::pair<int,int>> &edges) {
    Adj a{};
    for (int e = 0; e < int(edges.size()); ++e) if (graph & (Mask(1) << e)) {
        const auto [u,v] = edges[e]; a[u] |= Mask(1) << v; a[v] |= Mask(1) << u;
    }
    return a;
}

bool bipartite(const Adj &a, int n) {
    std::array<int,7> colors; colors.fill(-1);
    std::array<int,7> queue{};
    for (int root = 0; root < n; ++root) if (colors[root] == -1) {
        int head=0, tail=0; queue[tail++]=root; colors[root]=0;
        while (head < tail) {
            const int u=queue[head++];
            for (int v=0; v<n; ++v) if (a[u] & (Mask(1)<<v)) {
                if (colors[v] == colors[u]) return false;
                if (colors[v] == -1) {colors[v]=1-colors[u]; queue[tail++]=v;}
            }
        }
    }
    return true;
}

std::vector<std::vector<int>> components(const Adj &a, int n) {
    std::vector<std::vector<int>> result;
    Mask seen=0;
    for (int root=0; root<n; ++root) if (!(seen & (Mask(1)<<root))) {
        std::vector<int> part{root}; seen |= Mask(1)<<root;
        for (std::size_t j=0; j<part.size(); ++j)
            for(int v=0;v<n;++v)
                if ((a[part[j]] & (Mask(1)<<v)) && !(seen & (Mask(1)<<v))) {
                    seen |= Mask(1)<<v; part.push_back(v);
                }
        result.push_back(part);
    }
    return result;
}

bool affine_accept(const Adj &a, int n, bool omit_anchors=false,
                   bool single_global_factor=false, bool first_branch_only=false) {
    auto parts=components(a,n);
    if (single_global_factor && n) {
        std::vector<int> all(n); std::iota(all.begin(),all.end(),0); parts={all};
    }
    for (const auto &part : parts) {
        const int N=int(part.size()), d=2*N;
        std::vector<Row> equations;
        for(int i=0;i<N;++i) for(int j=i+1;j<N;++j) {
            const int u=part[i], v=part[j]; Row line=0;
            if(a[u] & (Mask(1)<<v)) line=(Row(1)<<i)|(Row(1)<<j)|(Row(1)<<d);
            for(int k=0;k<N;++k)
                if((a[u] & (Mask(1)<<part[k])) && (a[v] & (Mask(1)<<part[k])))
                    line |= Row(1)<<(N+k);
            equations.push_back(line);
        }
        Row solution=0;
        if (omit_anchors) {if(!solve(equations,d,solution)) return false; continue;}
        equations.push_back(Row(1)<<N); // u at the first vertex is zero.
        if(solve(equations,d,solution)) continue;
        if(first_branch_only) return false;
        Row anchor=0;
        for(int j=0;j<N;++j) if(a[part[0]] & (Mask(1)<<part[j])) anchor |= Row(1)<<(N+j);
        equations.back()=anchor; // (A u) at the first vertex is zero.
        if(!solve(equations,d,solution)) return false;
    }
    return true;
}

int main(int argc,char **argv) {
    if (argc != 3) {
        std::cerr << "Usage: graph_exhaustive MAX_ORDER OUTPUT.json\n";
        return 1;
    }
    const int max_n=std::stoi(argv[1]);
    if(max_n<0 || max_n>7) {
        std::cerr << "This finite audit supports max_n from 0 through 7.\n";
        return 1;
    }
    const auto started=std::chrono::steady_clock::now();
    const std::filesystem::path output(argv[2]);
    if (std::filesystem::exists(output)) {
        std::cerr << "Output already exists; select a fresh path.\n";
        return 1;
    }
    if (output.has_parent_path()) std::filesystem::create_directories(output.parent_path());
    std::ofstream report(output);
    if (!report) {
        std::cerr << "Cannot open output file.\n";
        return 1;
    }
    report << "{\n  \"method\": \"Independent complete labeled LC-orbit enumeration versus direct affine criterion\",\n  \"results\": [\n";
    for(int n=0;n<=max_n;++n) {
        std::vector<std::pair<int,int>> edges;
        for(int i=0;i<n;++i) for(int j=i+1;j<n;++j) edges.push_back({i,j});
        const Mask total=Mask(1)<<edges.size();
        std::vector<Mask> toggle(Mask(1)<<n,0);
        for(Mask set=0;set<(Mask(1)<<n);++set)
            for(int e=0;e<int(edges.size());++e)
                if((set & (Mask(1)<<edges[e].first)) && (set & (Mask(1)<<edges[e].second)))
                    toggle[set] |= Mask(1)<<e;
        std::vector<unsigned char> seen(total,0);
        std::uint64_t orbits=0, positives=0, negatives=0, mismatches=0;
        std::uint64_t omitted_failures=0, single_failures=0, first_failures=0;
        Mask omitted_example=0,single_example=0,first_example=0;
        for(Mask seed=0;seed<total;++seed) if(!seen[seed]) {
            ++orbits; std::vector<Mask> orbit{seed}; seen[seed]=1; bool has_bipartite=false;
            for(std::size_t p=0;p<orbit.size();++p) {
                const Mask graph=orbit[p]; const Adj a=adjacency(graph,edges);
                has_bipartite |= bipartite(a,n);
                for(int v=0;v<n;++v) {
                    const Mask next=graph ^ toggle[a[v]];
                    if(!seen[next]) {seen[next]=1; orbit.push_back(next);}
                }
            }
            if(has_bipartite) positives += orbit.size(); else negatives += orbit.size();
            for(Mask graph:orbit) {
                const Adj a=adjacency(graph,edges);
                if(affine_accept(a,n)!=has_bipartite) {
                    ++mismatches; std::cerr<<"MISMATCH n="<<n<<" graph="<<graph<<"\n";
                }
                if(n<=6) {
                    if(affine_accept(a,n,true)!=has_bipartite) {++omitted_failures;omitted_example=graph;}
                    if(affine_accept(a,n,false,true)!=has_bipartite) {++single_failures;single_example=graph;}
                    if(affine_accept(a,n,false,false,true)!=has_bipartite) {++first_failures;first_example=graph;}
                }
            }
        }
        std::cout<<"n="<<n<<" graphs="<<total<<" orbits="<<orbits<<" yes="<<positives
                 <<" no="<<negatives<<" mismatches="<<mismatches<<std::endl;
        report << (n?",\n":"") << "    {\"n\": "<<n<<", \"graphs\": "<<total
               <<", \"labeled_lc_orbits\": "<<orbits<<", \"yes\": "<<positives
               <<", \"no\": "<<negatives<<", \"mismatches\": "<<mismatches
               <<", \"mutations_checked\": "<<(n<=6?"true":"false")
               <<", \"without_anchors_failures\": "<<omitted_failures
               <<", \"without_anchors_last_graph\": "<<omitted_example
               <<", \"single_factor_failures\": "<<single_failures
               <<", \"single_factor_last_graph\": "<<single_example
               <<", \"first_branch_only_failures\": "<<first_failures
               <<", \"first_branch_only_last_graph\": "<<first_example<<"}";
        report.flush();
        if(mismatches) return 2;
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    report << "\n  ],\n  \"elapsed_seconds\": "<<seconds<<"\n}\n";
    std::cout<<"Elapsed seconds: "<<seconds<<std::endl;
}
