#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <time.h>
#include <ctime>
#include <string>
#include <string.h>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <vector>

using namespace std;

#define MAX_VAL 9999999
#define EPCION  0.001

#define EPSILON 0.9
#define GM 0.5
#define AF 0.5
#define FORGET_RATE 0.5

clock_t startTime;
double  bestTime;

int cuttingTime;
string inFile;
int p, r, verNum;

double gama = 0.09;
double alfa = 0.075;
double deta = 0.08;

double Sc = 0.0;
double Sb = MAX_VAL;

int far, far1;
double max_dis = 0.0;

int improve_flag_hubs = 0;

int W_param = 150;
int T_param = 5;
double q_ratio = 0.05;
double sl = 2.0;
double r0_param = 150.0;

double** edge;
double** temp_edge;
double** compute_edge;
double** compute_edge2;
double** compute_edge3;
int** seq_edge;
int** Path;
int** H;
int** cen_of_node;
int** node_of_cen;
int** idx_of_node;
int** idx_of_cen;
double** qtable;
double** reward;

double* edge_arr = nullptr;
double* temp_edge_arr = nullptr;
double* compute_edge_arr = nullptr;
double* compute_edge2_arr = nullptr;
double* compute_edge3_arr = nullptr;
double* qtable_arr = nullptr;
double* reward_arr = nullptr;
int* seq_edge_arr = nullptr;
int* Path_arr = nullptr;
int* H_arr = nullptr;
int* cen_of_node_arr = nullptr;
int* node_of_cen_arr = nullptr;
int* idx_of_node_arr = nullptr;
int* idx_of_cen_arr = nullptr;

int* cur_cen;
int* cur_idx;
int* cenFlag;
int* cen_len;
int* node_len;
int* U;
int* U1;
int* V;
int* V1;
int* qcenFlag;
int* hubs;
int* hubs_prev;
int* diff1;
int* diff2;

int criticalVer;
int criticalVer1;

int iter = 0;
int w_count = 0;

int random_int(int n) { return rand() % n; }

void freememo()
{
    delete[] edge_arr;
    delete[] temp_edge_arr;
    delete[] compute_edge_arr;
    delete[] compute_edge2_arr;
    delete[] compute_edge3_arr;
    delete[] qtable_arr;
    delete[] reward_arr;
    delete[] seq_edge_arr;
    delete[] Path_arr;
    delete[] H_arr;
    delete[] cen_of_node_arr;
    delete[] node_of_cen_arr;
    delete[] idx_of_node_arr;
    delete[] idx_of_cen_arr;

    delete[] edge;
    delete[] temp_edge;
    delete[] compute_edge;
    delete[] compute_edge2;
    delete[] compute_edge3;
    delete[] seq_edge;
    delete[] Path;
    delete[] H;
    delete[] cen_of_node;
    delete[] node_of_cen;
    delete[] idx_of_node;
    delete[] idx_of_cen;
    delete[] qtable;
    delete[] reward;

    delete[] cur_cen;
    delete[] cur_idx;
    delete[] cenFlag;
    delete[] cen_len;
    delete[] node_len;
    delete[] U;
    delete[] U1;
    delete[] V;
    delete[] V1;
    delete[] qcenFlag;
    delete[] hubs;
    delete[] hubs_prev;
    delete[] diff1;
    delete[] diff2;
}

void read_initial(const string& file)
{
    ifstream FIC(file);
    if (FIC.fail()) { exit(0); }

    FIC >> verNum;

    edge = new double* [verNum];
    temp_edge = new double* [verNum];
    compute_edge = new double* [verNum];
    compute_edge2 = new double* [verNum];
    compute_edge3 = new double* [verNum];
    seq_edge = new int* [verNum];
    Path = new int* [verNum];
    H = new int* [verNum];
    cen_of_node = new int* [verNum];
    node_of_cen = new int* [verNum];
    idx_of_node = new int* [verNum];
    idx_of_cen = new int* [verNum];
    qtable = new double* [verNum + 1];
    reward = new double* [verNum + 1];

    size_t n2 = (size_t)verNum * verNum;
    size_t nq = (size_t)(verNum + 1) * verNum;

    edge_arr = new double[n2]();
    temp_edge_arr = new double[n2]();
    compute_edge_arr = new double[n2]();
    compute_edge2_arr = new double[n2]();
    compute_edge3_arr = new double[n2]();
    qtable_arr = new double[nq]();
    reward_arr = new double[nq]();
    seq_edge_arr = new int[n2]();
    Path_arr = new int[n2]();
    H_arr = new int[n2]();
    cen_of_node_arr = new int[n2]();
    node_of_cen_arr = new int[n2]();
    idx_of_node_arr = new int[n2]();
    idx_of_cen_arr = new int[n2]();

    for (int i = 0; i < verNum; i++) {
        edge[i] = edge_arr + (size_t)i * verNum;
        temp_edge[i] = temp_edge_arr + (size_t)i * verNum;
        compute_edge[i] = compute_edge_arr + (size_t)i * verNum;
        compute_edge2[i] = compute_edge2_arr + (size_t)i * verNum;
        compute_edge3[i] = compute_edge3_arr + (size_t)i * verNum;
        seq_edge[i] = seq_edge_arr + (size_t)i * verNum;
        Path[i] = Path_arr + (size_t)i * verNum;
        H[i] = H_arr + (size_t)i * verNum;
        cen_of_node[i] = cen_of_node_arr + (size_t)i * verNum;
        node_of_cen[i] = node_of_cen_arr + (size_t)i * verNum;
        idx_of_node[i] = idx_of_node_arr + (size_t)i * verNum;
        idx_of_cen[i] = idx_of_cen_arr + (size_t)i * verNum;
        qtable[i] = qtable_arr + (size_t)i * verNum;
        reward[i] = reward_arr + (size_t)i * verNum;
    }
    qtable[verNum] = qtable_arr + (size_t)verNum * verNum;
    reward[verNum] = reward_arr + (size_t)verNum * verNum;

    cur_cen = new int[verNum];
    cur_idx = new int[verNum];
    cenFlag = new int[verNum];
    cen_len = new int[verNum];
    node_len = new int[verNum];
    U = new int[verNum];
    U1 = new int[verNum];
    V = new int[verNum];
    V1 = new int[verNum];
    qcenFlag = new int[verNum];
    hubs = new int[p];
    hubs_prev = new int[p];
    diff1 = new int[p];
    diff2 = new int[p];

    double* xCord = new double[verNum];
    double* yCord = new double[verNum];

    for (int i = 0; i < verNum; i++) {
        cenFlag[i] = 0;
        cen_len[i] = 0;
        node_len[i] = 0;
        xCord[i] = 0.0;
        yCord[i] = 0.0;
        for (int j = 0; j < verNum; j++) {
            edge[i][j] = 0.0;
            seq_edge[i][j] = j;
            H[i][j] = 0;
            Path[i][j] = -1;
            qtable[i][j] = 1.0;
            reward[i][j] = 1.0;
        }
    }
    for (int j = 0; j < verNum; j++) {
        qtable[verNum][j] = 1.0;
        reward[verNum][j] = 1.0;
    }

    string line;
    getline(FIC, line);
    getline(FIC, line);

    istringstream iss(line);
    vector<double> tokens;
    double tmp;
    while (iss >> tmp) tokens.push_back(tmp);
    int col_count = (int)tokens.size();

    FIC.close();
    FIC.open(file);
    FIC >> verNum;

    if (col_count == verNum) {
        for (int i = 0; i < verNum; i++)
            for (int j = 0; j < verNum; j++)
                FIC >> edge[i][j];
    }
    else if (col_count == 2) {
        for (int i = 0; i < verNum; i++)
            FIC >> xCord[i] >> yCord[i];
        for (int i = 0; i < verNum; i++)
            for (int j = 0; j < verNum; j++) {
                double dx = xCord[i] - xCord[j];
                double dy = yCord[i] - yCord[j];
                edge[i][j] = sqrt(dx * dx + dy * dy);
            }
    }
    else if (col_count == 3) {
        for (int i = 0; i < verNum; i++) {
            int idx_dummy;
            FIC >> idx_dummy >> xCord[i] >> yCord[i];
        }
        for (int i = 0; i < verNum; i++)
            for (int j = 0; j < verNum; j++) {
                double dx = xCord[i] - xCord[j];
                double dy = yCord[i] - yCord[j];
                edge[i][j] = sqrt(dx * dx + dy * dy);
            }
    }
    else {
        for (int i = 0; i < verNum; i++)
            for (int j = 0; j < verNum; j++)
                FIC >> edge[i][j];
    }

    FIC.close();

    for (int i = 0; i < verNum; i++)
        for (int j = 0; j < verNum; j++)
            temp_edge[i][j] = edge[i][j];

    for (int x = 0; x < verNum; x++)
        for (int i = 0; i < verNum - 1; i++)
            for (int j = 0; j < verNum - 1 - i; j++) {
                if (temp_edge[x][j] > temp_edge[x][j + 1]) {
                    double t1 = temp_edge[x][j];
                    temp_edge[x][j] = temp_edge[x][j + 1];
                    temp_edge[x][j + 1] = t1;
                    int t2 = seq_edge[x][j];
                    seq_edge[x][j] = seq_edge[x][j + 1];
                    seq_edge[x][j + 1] = t2;
                }
            }

    for (int i = 0; i < verNum; i++) {
        if (max_dis < temp_edge[i][verNum - 1]) {
            max_dis = temp_edge[i][verNum - 1];
            far = i;
            far1 = seq_edge[i][verNum - 1];
        }
    }

    delete[] xCord;
    delete[] yCord;
}

void write_outcome()
{
    ostringstream oss;
    oss << "p" << p << "r" << r << ".txt";
    ofstream outFile(oss.str(), ios_base::app);
    if (!outFile) { exit(0); }
    outFile << "p:" << p << " r:" << r
        << " cri:" << criticalVer << " cri1:" << criticalVer1
        << " time:" << bestTime << " Sb:" << Sb << endl;
    outFile.close();
}

bool restrict_nh()
{
    double temp_dis;
    double under1 = 0;
    double under2 = MAX_VAL;
    double under3 = MAX_VAL;

    if (!(cenFlag[criticalVer] || cenFlag[criticalVer1] || cenFlag[far] || cenFlag[far1])) {
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++) {
                double d2 = gama * edge[criticalVer][cur_cen[i]]
                    + alfa * edge[cur_cen[i]][cur_cen[j]]
                    + deta * edge[cur_cen[j]][criticalVer1];
                double d3 = gama * edge[far][cur_cen[i]]
                    + alfa * edge[cur_cen[i]][cur_cen[j]]
                    + deta * edge[cur_cen[j]][far1];
                if (d2 < under2) under2 = d2;
                if (d3 < under3) under3 = d3;
            }
    }
    else {
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
                compute_edge3[i][j] = alfa * edge[cur_cen[i]][cur_cen[j]];
        for (int k = 0; k < p; k++)
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++) {
                    double t = compute_edge3[i][k] + compute_edge3[k][j];
                    if (compute_edge3[i][j] > t) compute_edge3[i][j] = t;
                }
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
                if (compute_edge3[i][j] > under1) under1 = compute_edge3[i][j];

        if (cenFlag[criticalVer]) {
            if (cenFlag[criticalVer1])
                under2 = compute_edge3[cur_idx[criticalVer]][cur_idx[criticalVer1]];
            else
                for (int i = 0; i < p; i++) {
                    temp_dis = compute_edge3[cur_idx[criticalVer]][i]
                        + deta * edge[cur_cen[i]][criticalVer1];
                    if (temp_dis < under2) under2 = temp_dis;
                }
        }
        else {
            if (cenFlag[criticalVer1])
                for (int i = 0; i < p; i++) {
                    temp_dis = gama * edge[criticalVer][cur_cen[i]]
                        + compute_edge3[i][cur_idx[criticalVer1]];
                    if (temp_dis < under2) under2 = temp_dis;
                }
            else
                for (int i = 0; i < p; i++)
                    for (int j = 0; j < p; j++) {
                        temp_dis = gama * edge[criticalVer][cur_cen[i]]
                            + alfa * edge[cur_cen[i]][cur_cen[j]]
                            + deta * edge[cur_cen[j]][criticalVer1];
                        if (temp_dis < under2) under2 = temp_dis;
                    }
        }

        if (cenFlag[far]) {
            if (cenFlag[far1])
                under3 = compute_edge3[cur_idx[far]][cur_idx[far1]];
            else
                for (int i = 0; i < p; i++) {
                    temp_dis = compute_edge3[cur_idx[far]][i]
                        + deta * edge[cur_cen[i]][far1];
                    if (temp_dis < under3) under3 = temp_dis;
                }
        }
        else {
            if (cenFlag[far1])
                for (int i = 0; i < p; i++) {
                    temp_dis = gama * edge[far][cur_cen[i]]
                        + compute_edge3[i][cur_idx[far1]];
                    if (temp_dis < under3) under3 = temp_dis;
                }
            else
                for (int i = 0; i < p; i++)
                    for (int j = 0; j < p; j++) {
                        temp_dis = gama * edge[far][cur_cen[i]]
                            + alfa * edge[cur_cen[i]][cur_cen[j]]
                            + deta * edge[cur_cen[j]][far1];
                        if (temp_dis < under3) under3 = temp_dis;
                    }
        }
    }

    return (under1 < Sb&& under2 < Sb&& under3 < Sb);
}

double FLoyd()
{
    for (int i = 0; i < verNum; i++)
        for (int j = 0; j < verNum; j++)
            compute_edge[i][j] = compute_edge2[i][j] = MAX_VAL;

    for (int i = 0; i < p; i++)
        for (int j = 0; j < p; j++)
            compute_edge2[cur_cen[i]][cur_cen[j]] = alfa * edge[cur_cen[i]][cur_cen[j]];

    for (int k = 0; k < p; k++) {
        int tc = cur_cen[k];
        for (int j = 0; j < verNum; j++) {
            if (!cenFlag[j]) {
                for (int jr = 0; jr < r; jr++) {
                    int    tc2 = cen_of_node[j][jr];
                    double d = alfa * edge[tc][tc2] + deta * edge[tc2][j];
                    if (d < compute_edge[tc][j]) {
                        compute_edge[tc][j] = d;
                        Path[tc][j] = tc2;
                    }
                }
            }
            else {
                compute_edge[tc][j] = compute_edge[j][tc] = alfa * edge[tc][j];
                Path[tc][j] = j;
            }
        }
    }

    for (int i = 0; i < verNum; i++) {
        if (cenFlag[i]) continue;
        for (int k = 0; k < r; k++) {
            int tc = cen_of_node[i][k];
            for (int j = 0; j < verNum; j++) {
                double d = gama * edge[i][tc] + compute_edge[tc][j];
                if (d < compute_edge[i][j]) {
                    compute_edge[i][j] = d;
                    Path[i][j] = tc;
                }
            }
        }
    }

    for (int k = 0; k < p; k++)
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++) {
                double t = compute_edge2[cur_cen[i]][cur_cen[k]]
                    + compute_edge2[cur_cen[k]][cur_cen[j]];
                if (compute_edge2[cur_cen[i]][cur_cen[j]] > t)
                    compute_edge2[cur_cen[i]][cur_cen[j]] = t;
            }

    for (int i = 0; i < verNum; i++) {
        if (cenFlag[i]) continue;
        for (int k = 0; k < r; k++) {
            int tc = cen_of_node[i][k];
            for (int m = 0; m < p; m++) {
                int    tc2 = cur_cen[m];
                double d1 = gama * edge[i][tc] + compute_edge2[tc][tc2];
                if (d1 < compute_edge2[i][tc2]) compute_edge2[i][tc2] = d1;
                double d2 = deta * edge[tc][i] + compute_edge2[tc2][tc];
                if (compute_edge2[tc2][i] > d2) compute_edge2[tc2][i] = d2;
            }
        }
    }

    double temp_scc = 0;
    for (int i = 0; i < verNum; i++) {
        int ci = cenFlag[i];
        for (int j = 0; j < verNum; j++) {
            double dd = (ci | cenFlag[j]) ? compute_edge2[i][j] : compute_edge[i][j];
            if (dd > temp_scc) {
                temp_scc = dd;
                criticalVer = i;
                criticalVer1 = j;
            }
        }
    }
    return temp_scc;
}

void SWAP(int add, int drop)
{
    cur_cen[cur_idx[drop]] = add;
    cur_idx[add] = cur_idx[drop];
    cenFlag[drop] = 0;
    cenFlag[add] = 1;

    for (int i = 0; i < cen_len[drop]; i++) {
        int nd = node_of_cen[drop][i];
        node_of_cen[add][i] = nd;
        idx_of_cen[add][nd] = i;
        H[nd][add] = 1;
        H[nd][drop] = 0;
        int pos = idx_of_node[nd][drop];
        cen_of_node[nd][pos] = add;
        idx_of_node[nd][add] = pos;
    }
    cen_len[add] = cen_len[drop];
    cen_len[drop] = 0;
}

void INSERT(int nd, int org_cen, int tar_cen)
{
    H[nd][tar_cen] = 1;
    H[nd][org_cen] = 0;
    int pos = idx_of_node[nd][org_cen];
    cen_of_node[nd][pos] = tar_cen;
    idx_of_node[nd][tar_cen] = pos;

    node_of_cen[tar_cen][cen_len[tar_cen]] = nd;
    idx_of_cen[tar_cen][nd] = cen_len[tar_cen];
    cen_len[tar_cen]++;

    cen_len[org_cen]--;
    int p1 = idx_of_cen[org_cen][nd];
    int p2 = node_of_cen[org_cen][cen_len[org_cen]];
    node_of_cen[org_cen][p1] = p2;
    idx_of_cen[org_cen][p2] = p1;
}

double update_dis_insert(int cri, int org_cen, int tar_cen)
{
    int tc, tc2;
    double d;

    for (int i = 0; i < verNum; i++) {
        if (i == cri) continue;
        if (Path[cri][i] == org_cen) {
            double best = MAX_VAL;
            for (int j = 0; j < r; j++) {
                tc = cen_of_node[cri][j];
                d = gama * edge[cri][tc] + compute_edge[tc][i];
                if (d < best) {
                    best = d;
                    compute_edge[cri][i] = d;
                    Path[cri][i] = tc;
                }
            }
        }
        else {
            d = gama * edge[cri][tar_cen] + compute_edge[tar_cen][i];
            if (d < compute_edge[cri][i]) {
                compute_edge[cri][i] = d;
                Path[cri][i] = tar_cen;
            }
        }
    }

    for (int i = 0; i < verNum; i++)
        compute_edge[i][cri] = compute_edge2[i][cri] = compute_edge2[cri][i] = MAX_VAL;

    for (int k = 0; k < p; k++) {
        tc = cur_cen[k];
        for (int j = 0; j < r; j++) {
            tc2 = cen_of_node[cri][j];
            d = alfa * edge[tc][tc2] + deta * edge[tc2][cri];
            if (d < compute_edge[tc][cri]) {
                compute_edge[tc][cri] = d;
                Path[tc][cri] = tc2;
            }
        }
    }

    for (int i = 0; i < verNum; i++) {
        if (cenFlag[i]) continue;
        for (int k = 0; k < r; k++) {
            tc = cen_of_node[i][k];
            d = gama * edge[i][tc] + compute_edge[tc][cri];
            if (d < compute_edge[i][cri]) {
                compute_edge[i][cri] = d;
                Path[i][cri] = tc;
            }
        }
    }

    for (int j = 0; j < r; j++) {
        tc = cen_of_node[cri][j];
        for (int k = 0; k < p; k++) {
            tc2 = cur_cen[k];
            double d1 = gama * edge[cri][tc] + compute_edge2[tc][tc2];
            if (d1 < compute_edge2[cri][tc2]) compute_edge2[cri][tc2] = d1;
            double d2 = deta * edge[tc][cri] + compute_edge2[tc2][tc];
            if (compute_edge2[tc2][cri] > d2) compute_edge2[tc2][cri] = d2;
        }
    }

    double temp_scc = 0;
    for (int i = 0; i < verNum; i++) {
        int ci = cenFlag[i];
        for (int j = 0; j < verNum; j++) {
            double dd = (ci | cenFlag[j]) ? compute_edge2[i][j] : compute_edge[i][j];
            if (dd > temp_scc) {
                temp_scc = dd;
                criticalVer = i;
                criticalVer1 = j;
            }
        }
    }
    return temp_scc;
}

double estimate_dis(int cri, int cri1)
{
    double dis = MAX_VAL, temp_dis;
    for (int i = 0; i < r; i++)
        for (int j = 0; j < r; j++) {
            temp_dis = gama * edge[cri][cen_of_node[cri][i]]
                + alfa * edge[cen_of_node[cri][i]][cen_of_node[cri1][j]]
                + deta * edge[cri1][cen_of_node[cri1][j]];
            if (temp_dis < dis) dis = temp_dis;
        }
    return dis;
}

void update_best()
{
    if (Sc < Sb - EPCION) {
        clock_t endTime = clock();
        bestTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
        Sb = Sc;
    }
}

double compute_removal_cost(int node, int hub, int alloc_count)
{
    if (!H[node][hub]) return MAX_VAL;

    int original_pos = idx_of_node[node][hub];
    int last_hub = cen_of_node[node][alloc_count - 1];
    cen_of_node[node][original_pos] = last_hub;
    idx_of_node[node][last_hub] = original_pos;
    H[node][hub] = 0;

    int pos_in_cen = idx_of_cen[hub][node];
    int last_node = node_of_cen[hub][cen_len[hub] - 1];
    node_of_cen[hub][pos_in_cen] = last_node;
    idx_of_cen[hub][last_node] = pos_in_cen;
    cen_len[hub]--;

    double new_Sc = FLoyd();

    node_of_cen[hub][cen_len[hub]] = last_node;
    idx_of_cen[hub][last_node] = cen_len[hub];
    cen_len[hub]++;

    node_of_cen[hub][pos_in_cen] = node;
    idx_of_cen[hub][node] = pos_in_cen;

    cen_of_node[node][original_pos] = hub;
    idx_of_node[node][hub] = original_pos;
    H[node][hub] = 1;

    return new_Sc;
}

void greedy_removal()
{
    for (int node = 0; node < verNum; node++) {
        if (cenFlag[node]) continue;

        int alloc_count = 0;
        for (int i = 0; i < p; i++)
            if (H[node][cur_cen[i]]) alloc_count++;

        int to_remove = alloc_count - r;
        if (to_remove <= 0) continue;

        for (int del = 0; del < to_remove; del++) {
            int cur_alloc = alloc_count - del;
            double best_cost = MAX_VAL;
            int best_hub = -1;
            int best_pos = -1;

            for (int i = 0; i < cur_alloc; i++) {
                int    hub = cen_of_node[node][i];
                double cost = compute_removal_cost(node, hub, cur_alloc);
                if (cost < best_cost) {
                    best_cost = cost;
                    best_hub = hub;
                    best_pos = i;
                }
            }

            if (best_hub == -1) break;

            int last_hub = cen_of_node[node][cur_alloc - 1];
            cen_of_node[node][best_pos] = last_hub;
            idx_of_node[node][last_hub] = best_pos;
            H[node][best_hub] = 0;

            int pos_in_cen = idx_of_cen[best_hub][node];
            int last_node = node_of_cen[best_hub][cen_len[best_hub] - 1];
            node_of_cen[best_hub][pos_in_cen] = last_node;
            idx_of_cen[best_hub][last_node] = pos_in_cen;
            cen_len[best_hub]--;
        }
    }

    Sc = FLoyd();
}

bool na_descent()
{
    int org_cen, tar_cen, org_cen2, tar_cen2;
    int x, x1, y, y1, r0, r1, r20, r21;
    double delta;
    bool any_improve = false;

    int na_improve = 1;
    while (na_improve)
    {
        na_improve = 0;
        int cri = criticalVer;
        int cri1 = criticalVer1;

        if (!cenFlag[cri])
        {
            for (int i = 0; i < r; i++) U[i] = cen_of_node[cri][i];
            x = r;
            while (x > 0)
            {
                r0 = random_int(x);
                org_cen = U[r0];
                x1 = 0;
                for (int j = 0; j < p; j++)
                    if (!H[cri][cur_cen[j]]) U1[x1++] = cur_cen[j];

                while (x1 > 0)
                {
                    r1 = random_int(x1);
                    tar_cen = U1[r1];
                    if (gama * edge[cri][tar_cen] + compute_edge[tar_cen][cri1] < Sc)
                    {
                        INSERT(cri, org_cen, tar_cen);
                        delta = update_dis_insert(cri, org_cen, tar_cen) - Sc;
                        if (delta < -EPCION) {
                            Sc += delta;
                            any_improve = true;
                            na_improve = 1;
                            update_best();
                            goto next_vnd_round;
                        }
                        INSERT(cri, tar_cen, org_cen);
                        update_dis_insert(cri, tar_cen, org_cen);
                    }
                    x1--; U1[r1] = U1[x1];
                }
                x--; U[r0] = U[x];
            }
        }

        if (cri != cri1 && !cenFlag[cri1])
        {
            for (int i = 0; i < r; i++) U[i] = cen_of_node[cri1][i];
            x = r;
            while (x > 0)
            {
                r0 = random_int(x);
                org_cen = U[r0];
                x1 = 0;
                for (int j = 0; j < p; j++)
                    if (!H[cri1][cur_cen[j]]) U1[x1++] = cur_cen[j];

                while (x1 > 0)
                {
                    r1 = random_int(x1);
                    tar_cen = U1[r1];
                    if (compute_edge[cri][tar_cen] + deta * edge[tar_cen][cri1] < Sc)
                    {
                        INSERT(cri1, org_cen, tar_cen);
                        delta = update_dis_insert(cri1, org_cen, tar_cen) - Sc;
                        if (delta < -EPCION) {
                            Sc += delta;
                            any_improve = true;
                            na_improve = 1;
                            update_best();
                            goto next_vnd_round;
                        }
                        INSERT(cri1, tar_cen, org_cen);
                        update_dis_insert(cri1, tar_cen, org_cen);
                    }
                    x1--; U1[r1] = U1[x1];
                }
                x--; U[r0] = U[x];
            }
        }

        if (cri != cri1 && !cenFlag[cri] && !cenFlag[cri1])
        {
            for (int i = 0; i < r; i++) U[i] = cen_of_node[cri][i];
            x = r;
            while (x > 0)
            {
                r0 = random_int(x);
                org_cen = U[r0];

                for (int j = 0; j < r; j++) V[j] = cen_of_node[cri1][j];
                y = r;

                while (y > 0)
                {
                    r20 = random_int(y);
                    org_cen2 = V[r20];

                    int y1_base = 0;
                    for (int l = 0; l < p; l++)
                        if (!H[cri1][cur_cen[l]]) V1[y1_base++] = cur_cen[l];

                    x1 = 0;
                    for (int k = 0; k < p; k++)
                        if (!H[cri][cur_cen[k]]) U1[x1++] = cur_cen[k];

                    while (x1 > 0)
                    {
                        r1 = random_int(x1);
                        tar_cen = U1[r1];

                        INSERT(cri, org_cen, tar_cen);
                        bool cri_dis_updated = false;
                        y1 = y1_base;

                        while (y1 > 0)
                        {
                            r21 = random_int(y1);
                            tar_cen2 = V1[r21];

                            INSERT(cri1, org_cen2, tar_cen2);

                            if (estimate_dis(cri, cri1) < Sc)
                            {
                                if (!cri_dis_updated) {
                                    update_dis_insert(cri, org_cen, tar_cen);
                                    cri_dis_updated = true;
                                }
                                delta = update_dis_insert(cri1, tar_cen2, org_cen2) - Sc;
                                if (delta < -EPCION) {
                                    Sc += delta;
                                    any_improve = true;
                                    na_improve = 1;
                                    update_best();
                                    goto next_vnd_round;
                                }
                                INSERT(cri1, tar_cen2, org_cen2);
                                update_dis_insert(cri1, tar_cen2, org_cen2);
                            }
                            else {
                                INSERT(cri1, tar_cen2, org_cen2);
                            }

                            y1--; V1[r21] = V1[y1];
                        }

                        if (na_improve) break;

                        INSERT(cri, tar_cen, org_cen);
                        if (cri_dis_updated)
                            update_dis_insert(cri, tar_cen, org_cen);
                        x1--; U1[r1] = U1[x1];
                    }
                    if (na_improve) break;
                    y--; V[r20] = V[y];
                }
                if (na_improve) break;
                x--; U[r0] = U[x];
            }
        }

    next_vnd_round:;
    }

    return any_improve;
}

void initial_construct()
{
    for (int i = 0; i < verNum; i++) {
        cenFlag[i] = 0;
        cen_len[i] = 0;
        node_len[i] = 0;
        for (int j = 0; j < verNum; j++) H[i][j] = 0;
    }

    int clustNum = 0;
    while (clustNum < p) {
        int s = random_int(verNum);
        while (cenFlag[s]) s = random_int(verNum);
        cur_cen[clustNum] = s;
        cur_idx[s] = clustNum;
        cenFlag[s] = 1;
        clustNum++;
    }

    for (int i = 0; i < verNum; i++) {
        int u = random_int(p);
        int cnt = 0;
        for (int j = 0; j < p; j++) {
            int tc = cur_cen[(j + u) % p];
            if (tc == i) continue;
            node_of_cen[tc][cen_len[tc]] = i;
            idx_of_cen[tc][i] = cen_len[tc];
            cen_len[tc]++;
            cen_of_node[i][cnt] = tc;
            idx_of_node[i][tc] = cnt;
            H[i][tc] = 1;
            cnt++;
            if (cnt == r) break;
        }
    }
}

void reward_updating()
{
    for (int i = 0; i < p; i++) {
        if (hubs_prev[i] < 0 || hubs_prev[i] >= verNum) return;
        if (cur_cen[i] < 0 || cur_cen[i] >= verNum) return;
    }

    for (int i = 0; i < p; i++) {
        bool unchanged = false;
        for (int j = 0; j < p; j++)
            if (hubs_prev[i] == cur_cen[j]) { unchanged = true; break; }

        if (unchanged) {
            for (int s = 0; s <= verNum; s++) {
                reward[s][hubs_prev[i]] += sl;
                if (reward[s][hubs_prev[i]] > r0_param)
                    reward[s][hubs_prev[i]] *= FORGET_RATE;
            }
        }
        else {
            for (int s = 0; s <= verNum; s++) {
                reward[s][hubs_prev[i]] -= sl;
                if (reward[s][hubs_prev[i]] < 0.1)
                    reward[s][hubs_prev[i]] = 0.1;
            }
        }
    }

    for (int i = 0; i < p; i++) {
        bool in_prev = false;
        for (int j = 0; j < p; j++)
            if (cur_cen[i] == hubs_prev[j]) { in_prev = true; break; }
        if (!in_prev) {
            for (int s = 0; s <= verNum; s++) {
                reward[s][cur_cen[i]] += sl;
                if (reward[s][cur_cen[i]] > r0_param)
                    reward[s][cur_cen[i]] *= FORGET_RATE;
            }
        }
    }

    if (Sc < Sb - EPCION) {
        for (int i = 0; i < p; i++)
            for (int s = 0; s <= verNum; s++) {
                reward[s][cur_cen[i]] += sl;
                if (reward[s][cur_cen[i]] > r0_param)
                    reward[s][cur_cen[i]] *= FORGET_RATE;
            }
    }
}

void ql_construct()
{
    int    st = 0;
    int    st_prev = 0;
    double r_num = (double)rand() / RAND_MAX;

    for (int i = 0; i < p; i++) qcenFlag[cur_cen[i]] = 0;

    for (int i = 0; i < p; i++) {
        int    select_hub = -1;
        double max_q = -1.0;

        if (r_num < EPSILON) {
            for (int j = 0; j < verNum; j++)
                if (!qcenFlag[j] && qtable[st][j] > max_q) {
                    max_q = qtable[st][j];
                    select_hub = j;
                }
        }
        else {
            select_hub = random_int(verNum);
            while (qcenFlag[select_hub]) select_hub = random_int(verNum);
        }

        if (select_hub < 0) {
            select_hub = random_int(verNum);
            while (qcenFlag[select_hub]) select_hub = random_int(verNum);
        }

        st_prev = st;
        st = select_hub + 1;
        qcenFlag[select_hub] = 1;
        hubs[i] = select_hub;

        double next_max_q = 0;
        for (int j = 0; j < verNum; j++)
            if (qtable[st][j] > next_max_q) next_max_q = qtable[st][j];

        qtable[st_prev][select_hub] =
            (1.0 - GM) * qtable[st_prev][select_hub]
            + AF * (reward[st_prev][select_hub] + GM * next_max_q);
    }

    for (int i = 0; i < p; i++) qcenFlag[hubs[i]] = 0;

    int diff_num = 0, drop_num = 0;
    for (int i = 0; i < p; i++) {
        bool found = false;
        for (int j = 0; j < p; j++)
            if (hubs[i] == cur_cen[j]) { found = true; break; }
        if (!found) diff1[diff_num++] = hubs[i];
    }
    for (int i = 0; i < p; i++) {
        bool found = false;
        for (int j = 0; j < p; j++)
            if (cur_cen[i] == hubs[j]) { found = true; break; }
        if (!found) diff2[drop_num++] = cur_cen[i];
    }
    for (int i = 0; i < diff_num; i++)
        SWAP(diff1[i], diff2[i]);

    greedy_removal();

    Sc = FLoyd();
    for (int i = 0; i < p; i++) hubs_prev[i] = cur_cen[i];
    update_best();
}

void threshold_search()
{
    int    temp_id, idx;
    double delta;

    for (int t = 0; t < T_param; t++)
    {
        bool nh_accepted = false;

        for (int i = 0; i < verNum; i++) U[i] = i;
        int x = verNum;

        while (x > 0 && !nh_accepted)
        {
            int r0 = random_int(x);
            temp_id = U[r0];

            if (cenFlag[temp_id]) { x--; U[r0] = U[x]; continue; }

            for (int j = 0; j < p; j++) U1[j] = cur_cen[j];
            int x1 = p;

            while (x1 > 0 && !nh_accepted)
            {
                int r1 = random_int(x1);
                idx = U1[r1];

                int saved_cri = criticalVer;
                int saved_cri1 = criticalVer1;

                SWAP(temp_id, idx);

                if (restrict_nh()) {
                    delta = FLoyd() - Sc;

                    if (delta < -EPCION) {
                        Sc += delta;
                        improve_flag_hubs = 1;
                        nh_accepted = true;
                        update_best();
                    }
                    else if (Sc + delta <= (1.0 + q_ratio) * Sb) {
                        Sc += delta;
                        improve_flag_hubs = 0;
                        nh_accepted = true;
                    }
                    else {
                        SWAP(idx, temp_id);
                        criticalVer = saved_cri;
                        criticalVer1 = saved_cri1;
                    }
                }
                else {
                    SWAP(idx, temp_id);
                    criticalVer = saved_cri;
                    criticalVer1 = saved_cri1;
                }

                if (!nh_accepted) { x1--; U1[r1] = U1[x1]; }
            }
            if (!nh_accepted) { x--; U[r0] = U[x]; }
        }

        if (!nh_accepted) {
            reward_updating();
            ql_construct();
            na_descent();
            w_count = 0;
            return;
        }

        na_descent();
    }

    if (Sc < Sb - EPCION)
        w_count = 0;
    else
        w_count++;

    if (w_count >= W_param) {
        reward_updating();
        ql_construct();
        na_descent();
        w_count = 0;
    }
}

int main(int argc, char* argv[])
{
    if (argc > 1) inFile = argv[1];
    if (argc > 2) p = atoi(argv[2]);
    if (argc > 3) r = atoi(argv[3]);
    if (argc > 4) cuttingTime = atoi(argv[4]);
    if (argc > 5) gama = atof(argv[5]);
    if (argc > 6) alfa = atof(argv[6]);
    if (argc > 7) deta = atof(argv[7]);

    srand((unsigned)time(NULL));
    startTime = clock();

    read_initial(inFile);

    initial_construct();
    Sc = FLoyd();
    Sb = Sc;
    for (int i = 0; i < p; i++) hubs_prev[i] = cur_cen[i];

    na_descent();
    Sb = Sc;

    while (true)
    {
        double elapsed = (double)(clock() - startTime) / CLOCKS_PER_SEC;
        if (elapsed >= cuttingTime) break;
        iter++;
        threshold_search();
    }

    write_outcome();

    freememo();
    return 0;
}
