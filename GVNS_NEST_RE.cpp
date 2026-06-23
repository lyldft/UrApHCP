#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iomanip>

using namespace std;

#define MAX_VAL 9999999
#define EPCION  0.001

clock_t startTime;
int cuttingTime;
string inFile;
int p;
int r;

double gama = 1;
double alfa = 0.75;
double deta = 1;

int verNum;
double** edge;
double** compute_edge;
double** compute_edge2;
int** H;
int* cur_cen;
int* cenFlag;
int* cur_idx;
int** cen_of_node;
int** idx_of_node;
int** node_of_cen;
int** idx_of_cen;
int* cen_len;
int* U;

double Sc = 0.0;
double Sb = MAX_VAL;

inline bool time_up() {
    return ((double)(clock() - startTime) / CLOCKS_PER_SEC >= cuttingTime);
}

int random_int(int n) {
    if (n <= 0) return 0;
    return rand() % n;
}

void read_initial(string fileName) {
    ifstream FIC(fileName);
    if (!FIC) {
        cout << "### Error open, File_Name " << fileName << endl;
        exit(0);
    }

    FIC >> verNum;

    edge = new double* [verNum];
    compute_edge = new double* [verNum];
    compute_edge2 = new double* [verNum];
    H = new int* [verNum];
    cur_cen = new int[verNum];
    cenFlag = new int[verNum];
    cur_idx = new int[verNum];
    cen_of_node = new int* [verNum];
    idx_of_node = new int* [verNum];
    node_of_cen = new int* [verNum];
    idx_of_cen = new int* [verNum];
    cen_len = new int[verNum];
    U = new int[verNum];

    for (int i = 0; i < verNum; i++) {
        edge[i] = new double[verNum];
        compute_edge[i] = new double[verNum];
        compute_edge2[i] = new double[verNum];
        H[i] = new int[verNum];
        cen_of_node[i] = new int[verNum];
        idx_of_node[i] = new int[verNum];
        node_of_cen[i] = new int[verNum];
        idx_of_cen[i] = new int[verNum];

        cur_cen[i] = 0;
        cenFlag[i] = 0;
        cur_idx[i] = 0;
        cen_len[i] = 0;

        for (int j = 0; j < verNum; j++) {
            edge[i][j] = 0.0;
            compute_edge[i][j] = 0.0;
            compute_edge2[i][j] = 0.0;
            H[i][j] = 0;
            cen_of_node[i][j] = 0;
            idx_of_node[i][j] = 0;
            node_of_cen[i][j] = 0;
            idx_of_cen[i][j] = 0;
        }
    }

    string firstLine;
    getline(FIC, firstLine);
    while (firstLine.empty() && getline(FIC, firstLine));

    istringstream probe(firstLine);
    vector<double> tokens;
    double tv;
    while (probe >> tv) tokens.push_back(tv);

    bool coordMode = (tokens.size() == 2);
    double* xCoord = new double[verNum];
    double* yCoord = new double[verNum];

    if (coordMode) {
        xCoord[0] = tokens[0];
        yCoord[0] = tokens[1];
        for (int i = 1; i < verNum; i++) FIC >> xCoord[i] >> yCoord[i];
        for (int i = 0; i < verNum; i++)
            for (int j = 0; j < verNum; j++) {
                if (i == j) { edge[i][j] = 0.0; continue; }
                double dx = xCoord[i] - xCoord[j];
                double dy = yCoord[i] - yCoord[j];
                edge[i][j] = sqrt(dx * dx + dy * dy);
            }
    }
    else {
        for (int j = 0; j < (int)tokens.size(); j++) edge[0][j] = tokens[j];
        for (int j = (int)tokens.size(); j < verNum; j++) FIC >> edge[0][j];
        for (int i = 1; i < verNum; i++)
            for (int j = 0; j < verNum; j++) FIC >> edge[i][j];
    }

    FIC.close();
    delete[] xCoord;
    delete[] yCoord;
}

void save_hubs(int* buf) {
    for (int i = 0; i < p; i++) buf[i] = cur_cen[i];
}

void rebuild_flags() {
    for (int i = 0; i < verNum; i++) { cenFlag[i] = 0; cur_idx[i] = 0; }
    for (int i = 0; i < p; i++) { cenFlag[cur_cen[i]] = 1; cur_idx[cur_cen[i]] = i; }
}

void restore_hubs(int* buf) {
    for (int i = 0; i < p; i++) cur_cen[i] = buf[i];
    rebuild_flags();
}

void clear_assignments() {
    for (int i = 0; i < verNum; i++) {
        cen_len[i] = 0;
        for (int j = 0; j < verNum; j++) {
            H[i][j] = 0;
            cen_of_node[i][j] = 0;
            idx_of_node[i][j] = 0;
            node_of_cen[i][j] = 0;
            idx_of_cen[i][j] = 0;
        }
    }
}

void SWAP(int add, int drop) {
    int pos = cur_idx[drop];
    cur_cen[pos] = add;
    cur_idx[add] = pos;
    cenFlag[drop] = 0;
    cenFlag[add] = 1;
}

void INSERT(int nd, int org_cen, int tar_cen) {
    H[nd][tar_cen] = 1;
    H[nd][org_cen] = 0;

    int pos = idx_of_node[nd][org_cen];
    cen_of_node[nd][pos] = tar_cen;
    idx_of_node[nd][tar_cen] = pos;

    node_of_cen[tar_cen][cen_len[tar_cen]] = nd;
    idx_of_cen[tar_cen][nd] = cen_len[tar_cen];
    cen_len[tar_cen]++;

    int pos2 = idx_of_cen[org_cen][nd];
    int last = node_of_cen[org_cen][cen_len[org_cen] - 1];
    node_of_cen[org_cen][pos2] = last;
    idx_of_cen[org_cen][last] = pos2;
    cen_len[org_cen]--;
}

void realloc_random() {
    clear_assignments();

    for (int nd = 0; nd < verNum; nd++) {
        if (time_up()) return;
        if (!cenFlag[nd]) continue;

        int cnt = 0;
        cen_of_node[nd][cnt] = nd;
        idx_of_node[nd][nd] = cnt;
        H[nd][nd] = 1;
        node_of_cen[nd][cen_len[nd]] = nd;
        idx_of_cen[nd][nd] = cen_len[nd];
        cen_len[nd]++;
        cnt++;

        int guard = 0;
        while (cnt < r && guard < verNum * 4 && !time_up()) {
            int h = cur_cen[random_int(p)];
            guard++;
            if (h == nd || H[nd][h]) continue;
            cen_of_node[nd][cnt] = h; idx_of_node[nd][h] = cnt; H[nd][h] = 1;
            node_of_cen[h][cen_len[h]] = nd; idx_of_cen[h][nd] = cen_len[h]; cen_len[h]++;
            cnt++;
        }
        for (int k = 0; cnt < r && k < p; k++) {
            int h = cur_cen[k];
            if (h == nd || H[nd][h]) continue;
            cen_of_node[nd][cnt] = h; idx_of_node[nd][h] = cnt; H[nd][h] = 1;
            node_of_cen[h][cen_len[h]] = nd; idx_of_cen[h][nd] = cen_len[h]; cen_len[h]++;
            cnt++;
        }
    }

    for (int nd = 0; nd < verNum; nd++) {
        if (time_up()) return;
        if (cenFlag[nd]) continue;

        int cnt = 0, guard = 0;
        while (cnt < r && guard < verNum * 4 && !time_up()) {
            int h = cur_cen[random_int(p)];
            guard++;
            if (H[nd][h]) continue;
            cen_of_node[nd][cnt] = h; idx_of_node[nd][h] = cnt; H[nd][h] = 1;
            node_of_cen[h][cen_len[h]] = nd; idx_of_cen[h][nd] = cen_len[h]; cen_len[h]++;
            cnt++;
        }
        for (int k = 0; cnt < r && k < p; k++) {
            int h = cur_cen[k];
            if (H[nd][h]) continue;
            cen_of_node[nd][cnt] = h; idx_of_node[nd][h] = cnt; H[nd][h] = 1;
            node_of_cen[h][cen_len[h]] = nd; idx_of_cen[h][nd] = cen_len[h]; cen_len[h]++;
            cnt++;
        }
    }
}

void realloc_greedy() {
    clear_assignments();

    for (int nd = 0; nd < verNum; nd++) {
        if (time_up()) return;
        if (!cenFlag[nd]) continue;

        int cnt = 0;
        cen_of_node[nd][cnt] = nd; idx_of_node[nd][nd] = cnt; H[nd][nd] = 1;
        node_of_cen[nd][cen_len[nd]] = nd; idx_of_cen[nd][nd] = cen_len[nd]; cen_len[nd]++;
        cnt++;

        for (int k = 0; k < p && cnt < r; k++) {
            int h = cur_cen[k];
            if (h == nd || H[nd][h]) continue;
            cen_of_node[nd][cnt] = h; idx_of_node[nd][h] = cnt; H[nd][h] = 1;
            node_of_cen[h][cen_len[h]] = nd; idx_of_cen[h][nd] = cen_len[h]; cen_len[h]++;
            cnt++;
        }
    }

    for (int nd = 0; nd < verNum; nd++) {
        if (time_up()) return;
        if (cenFlag[nd]) continue;

        double* aloc_val = new double[p];
        int* order = new int[p];

        for (int k = 0; k < p; k++) {
            int h = cur_cen[k];
            double max_hh = 0.0;
            for (int k2 = 0; k2 < p; k2++) {
                double d = alfa * edge[h][cur_cen[k2]];
                if (d > max_hh) max_hh = d;
            }
            aloc_val[k] = gama * edge[nd][h] + max_hh;
            order[k] = k;
        }

        for (int i = 0; i < r; i++) {
            for (int j = i + 1; j < p; j++)
                if (aloc_val[order[j]] < aloc_val[order[i]]) swap(order[i], order[j]);
            int h = cur_cen[order[i]];
            cen_of_node[nd][i] = h; idx_of_node[nd][h] = i; H[nd][h] = 1;
            node_of_cen[h][cen_len[h]] = nd; idx_of_cen[h][nd] = cen_len[h]; cen_len[h]++;
        }

        delete[] aloc_val;
        delete[] order;
    }
}

double FLoyd() {
    for (int i = 0; i < verNum; i++) {
        if (time_up()) return Sc;
        for (int j = 0; j < verNum; j++) {
            compute_edge[i][j] = MAX_VAL;
            compute_edge2[i][j] = MAX_VAL;
        }
    }

    for (int i = 0; i < p; i++)
        for (int j = 0; j < p; j++)
            compute_edge2[cur_cen[i]][cur_cen[j]] = alfa * edge[cur_cen[i]][cur_cen[j]];

    for (int k = 0; k < p; k++) {
        if (time_up()) return Sc;
        int tc = cur_cen[k];
        for (int j = 0; j < verNum; j++) {
            if (!cenFlag[j]) {
                for (int jr = 0; jr < r; jr++) {
                    int tc2 = cen_of_node[j][jr];
                    double d = alfa * edge[tc][tc2] + deta * edge[tc2][j];
                    if (d < compute_edge[tc][j]) compute_edge[tc][j] = d;
                }
            }
            else {
                compute_edge[tc][j] = alfa * edge[tc][j];
            }
        }
    }

    for (int i = 0; i < verNum; i++) {
        if (time_up()) return Sc;
        if (cenFlag[i]) continue;
        for (int k = 0; k < r; k++) {
            int tc = cen_of_node[i][k];
            for (int j = 0; j < verNum; j++) {
                double d = gama * edge[i][tc] + compute_edge[tc][j];
                if (d < compute_edge[i][j]) compute_edge[i][j] = d;
            }
        }
    }

    for (int k = 0; k < p; k++) {
        if (time_up()) return Sc;
        int hk = cur_cen[k];
        for (int i = 0; i < p; i++) {
            int hi = cur_cen[i];
            for (int j = 0; j < p; j++) {
                int hj = cur_cen[j];
                double tx = compute_edge2[hi][hk] + compute_edge2[hk][hj];
                if (compute_edge2[hi][hj] > tx) compute_edge2[hi][hj] = tx;
            }
        }
    }

    for (int i = 0; i < verNum; i++) {
        if (time_up()) return Sc;
        if (cenFlag[i]) continue;
        for (int j = 0; j < r; j++) {
            int tc = cen_of_node[i][j];
            for (int k = 0; k < p; k++) {
                int tc2 = cur_cen[k];
                double d1 = gama * edge[i][tc] + compute_edge2[tc][tc2];
                if (d1 < compute_edge2[i][tc2]) compute_edge2[i][tc2] = d1;
                double d2 = deta * edge[tc][i] + compute_edge2[tc2][tc];
                if (compute_edge2[tc2][i] > d2) compute_edge2[tc2][i] = d2;
            }
        }
    }

    double temp_scc = 0.0;
    for (int i = 0; i < verNum; i++) {
        if (time_up()) return Sc;
        for (int j = 0; j < verNum; j++) {
            double val = (cenFlag[i] || cenFlag[j]) ? compute_edge2[i][j] : compute_edge[i][j];
            if (val > temp_scc) temp_scc = val;
        }
    }

    return temp_scc;
}

void init_construct() {
    for (int i = 0; i < verNum; i++) { cenFlag[i] = 0; cur_idx[i] = 0; cur_cen[i] = 0; }

    int selected = 0, guard = 0;
    while (selected < p && guard < verNum * 20 && !time_up()) {
        int h = random_int(verNum); guard++;
        if (cenFlag[h]) continue;
        cur_cen[selected] = h; cur_idx[h] = selected; cenFlag[h] = 1; selected++;
    }
    for (int i = 0; selected < p && i < verNum; i++) {
        if (!cenFlag[i]) { cur_cen[selected] = i; cur_idx[i] = selected; cenFlag[i] = 1; selected++; }
    }

    realloc_random();
}

bool na_first_improve() {
    for (int nd = 0; nd < verNum; nd++) {
        if (time_up()) return false;
        if (cenFlag[nd]) continue;
        for (int ri = 0; ri < r; ri++) {
            if (time_up()) return false;
            int org_cen = cen_of_node[nd][ri];
            for (int k = 0; k < p; k++) {
                if (time_up()) return false;
                int tar_cen = cur_cen[k];
                if (tar_cen == org_cen || H[nd][tar_cen]) continue;

                INSERT(nd, org_cen, tar_cen);
                double new_val = FLoyd();
                if (time_up()) { INSERT(nd, tar_cen, org_cen); return false; }
                if (new_val < Sc - EPCION) { Sc = new_val; return true; }
                else INSERT(nd, tar_cen, org_cen);
            }
        }
    }
    return false;
}

void na_local_descent() {
    bool improved = true;
    while (improved && !time_up()) improved = na_first_improve();
}

void restore_full(int* saved_hubs, double saved_Sc) {
    restore_hubs(saved_hubs);
    realloc_greedy();
    Sc = saved_Sc;
}

bool nestVND_once() {
    int* saved_hubs = new int[p];

    for (int add_hub = 0; add_hub < verNum; add_hub++) {
        if (time_up()) { delete[] saved_hubs; return false; }
        if (cenFlag[add_hub]) continue;

        for (int j = 0; j < p; j++) {
            if (time_up()) { delete[] saved_hubs; return false; }

            int drop_hub = cur_cen[j];
            save_hubs(saved_hubs);
            double saved_Sc = Sc;

            SWAP(add_hub, drop_hub);
            realloc_greedy();
            if (time_up()) { restore_full(saved_hubs, saved_Sc); delete[] saved_hubs; return false; }

            Sc = FLoyd();
            if (time_up()) { restore_full(saved_hubs, saved_Sc); delete[] saved_hubs; return false; }

            na_local_descent();
            if (time_up()) { restore_full(saved_hubs, saved_Sc); delete[] saved_hubs; return false; }

            double new_val = FLoyd();
            if (time_up()) { restore_full(saved_hubs, saved_Sc); delete[] saved_hubs; return false; }

            if (new_val < saved_Sc - EPCION) {
                Sc = new_val;
                if (Sc < Sb - EPCION) Sb = Sc;
                delete[] saved_hubs;
                return true;
            }
            else {
                restore_full(saved_hubs, saved_Sc);
            }
        }
    }

    delete[] saved_hubs;
    return false;
}

void nested_VND_loop() {
    bool improved = true;
    while (improved && !time_up()) improved = nestVND_once();
}

void shake(int k) {
    for (int i = 0; i < k; i++) {
        if (time_up()) return;
        int add_hub = random_int(verNum);
        int guard = 0;
        while (cenFlag[add_hub] && guard < verNum * 4 && !time_up()) {
            add_hub = random_int(verNum); guard++;
        }
        if (cenFlag[add_hub]) return;
        SWAP(add_hub, cur_cen[random_int(p)]);
    }
    realloc_random();
    if (!time_up()) Sc = FLoyd();
}

void gvns_nest() {
    init_construct();
    Sc = FLoyd();
    Sb = Sc;

    int kmax = 1;
    int* saved_hubs = new int[p];

    while (!time_up()) {
        nested_VND_loop();
        if (time_up()) break;

        if (Sc < Sb - EPCION) Sb = Sc;

        save_hubs(saved_hubs);
        double saved_Sc = Sc;

        int k = 1;
        while (k <= kmax && !time_up()) {
            restore_full(saved_hubs, saved_Sc);
            if (time_up()) break;
            shake(k);
            if (time_up()) break;
            nested_VND_loop();
            if (time_up()) break;

            if (Sc < saved_Sc - EPCION) {
                save_hubs(saved_hubs);
                saved_Sc = Sc;
                if (Sc < Sb - EPCION) Sb = Sc;
                k = 1;
            }
            else {
                k++;
            }
        }

        restore_full(saved_hubs, saved_Sc);
    }

    delete[] saved_hubs;
}

void write_outcome() {
    ostringstream outfilename;
    outfilename << "p" << p << "r" << r << "_GVNS_NEST_RANDOM.txt";
    ofstream outFile(outfilename.str(), ios_base::app);
    if (!outFile) { cerr << "无法打开文件" << endl; exit(0); }
    outFile << "p:" << p << " r:" << r << " Sb:" << fixed << setprecision(2) << Sb << endl;
    outFile.close();
}

void freememo() {
    for (int i = 0; i < verNum; i++) {
        delete[] edge[i];        delete[] compute_edge[i];
        delete[] compute_edge2[i]; delete[] H[i];
        delete[] cen_of_node[i]; delete[] idx_of_node[i];
        delete[] node_of_cen[i]; delete[] idx_of_cen[i];
    }
    delete[] edge;         delete[] compute_edge;
    delete[] compute_edge2; delete[] H;
    delete[] cur_cen;      delete[] cenFlag;
    delete[] cur_idx;      delete[] cen_of_node;
    delete[] idx_of_node;  delete[] node_of_cen;
    delete[] idx_of_cen;   delete[] cen_len;
    delete[] U;
}

int main(int argc, char** argv) {
    if (argc > 1) inFile = argv[1];
    if (argc > 2) p = atoi(argv[2]);
    if (argc > 3) r = atoi(argv[3]);
    if (argc > 4) cuttingTime = atoi(argv[4]);
    if (argc > 5) gama = atof(argv[5]);
    if (argc > 6) alfa = atof(argv[6]);
    if (argc > 7) deta = atof(argv[7]);

    cout << fixed << setprecision(2);
    srand((unsigned)time(NULL));

    read_initial(inFile);

    cout << "File:" << inFile << "  p=" << p << "  r=" << r
        << "  gama=" << gama << "  alfa=" << alfa << "  deta=" << deta << endl;

    startTime = clock();
    gvns_nest();

    cout << "Sb: " << Sb << endl;

    write_outcome();
    freememo();

    return 0;
}
