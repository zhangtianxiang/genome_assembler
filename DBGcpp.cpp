#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace std;

FILE *f_g, *f_s, *f_r;

int N, maxd = 50000;
vector<vector<int> > g;
vector<bool> visited, visited_back;
vector<int> maxto, maxto_back;
vector<int> S;

int dfs(int u, int dep = 0) {
    visited[u] = true;
    int maxdeepu = 1;
    if (dep == maxd) return maxdeepu;
    for (int v : g[u]) {
        if (visited[v]) continue;
        int maxdeepv = dfs(v, dep + 1);
        if (maxdeepv + 1 > maxdeepu) {
            maxdeepu = maxdeepv + 1;
            maxto[u] = v;
        }
    }
    return maxdeepu;
}

int main() {
    // open file
    f_g = fopen("Graph.txt", "r");
    assert(f_g != NULL);
    f_s = fopen("Start.txt", "r");
    assert(f_s != NULL);
    f_r = fopen("Result.txt", "w");
    assert(f_r != NULL);
    // read data
    fscanf(f_g, "%d", &N);
    int tote = 0;
    for (int i = 0; i < N; i++) {
        vector<int> now;
        int totnow, v;
        fscanf(f_g, "%d", &totnow);
        tote += totnow;
        for (int j = 0; j < totnow; j++) {
            fscanf(f_g, "%d", &v);
            now.push_back(v);
        }
        g.push_back(now);
        maxto_back.push_back(-1);
        visited_back.push_back(false);
    }
    fclose(f_g);
    // read st
    int tots;
    fscanf(f_s, "%d", &tots);
    for (int i = 0, j; i < tots; i++) {
        fscanf(f_s, "%d", &j);
        S.push_back(j);
    }
    fclose(f_s);
    // try dfs
    for (int s : S) {
        visited = visited_back;
        maxto = maxto_back;
        dfs(s);
        s = maxto[s];
        while (s != -1) {
            fprintf(f_r, "%d ", s);
            s = maxto[s];
        }
        fprintf(f_r, "\n");
    }
    fclose(f_r);
    return 0;
}