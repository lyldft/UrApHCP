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

int random_int(int n) { return rand() % n; }

void read_initial(string fileName) {
    ifstream FIC(fileName);
    if (!FIC) {
        cout << "### Erreur open, File_Name " << fileName << endl;
        exit(0);
    }

    FIC >> verNum;

    edge         = new double*[verNum];
    compute_edge = new double*[verNum];
    H            = new int*[verNum];
    cur_cen      = new int[verNum];
    cenFlag      = new int[verNum];
    cur_idx      = new int[verNum];
    cen_of_node  = new int*[verNum];
    idx_of_node  = new int*[verNum];
    node_of_cen  = new int*[verNum];
    idx_of_cen   = new int*[verNum];
    cen_len      = new int[verNum];
    U            = new int[verNum];

    for (int i = 0; i < verNum; i++) {
        edge[i]        = new double[verNum];
        compute_edge[i]= new double[verNum];
        H[i]           = new int[verNum];
        cen_of_node[i] = new int[verNum];
        idx_of_node[i] = new int[verNum];
        node_of_cen[i] = new int[verNum];
        idx_of_cen[i]  = new int[verNum];

        cenFlag[i] = 0;
        cen_len[i] = 0;

        for (int j = 0; j < verNum; j++) {
            edge[i][j]        = 0;
            compute_edge[i][j]= 0;
            H[i][j]           = 0;
        }
    }

    string firstLine;
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
            for (int j = 0; j < verNum; j++)
                edge[i][j] = (i == j) ? 0.0 :
                    sqrt((xCoord[i]-xCoord[j])*(xCoord[i]-xCoord[j]) +
                         (yCoord[i]-yCoord[j])*(yCoord[i]-yCoord[j]));
    } else {
        for (int j = 0; j < (int)tokens.size(); j++)
            edge[0][j] = tokens[j];
        for (int j = (int)tokens.size(); j < verNum; j++) FIC >> edge[0][j];
        for (int i = 1; i < verNum; i++)
            for (int j = 0; j < verNum; j++) FIC >> edge[i][j];
    }

    FIC.close();
    delete[] xCoord;
    delete[] yCoord;
}

void SWAP(int add, int drop) {
    cur_cen[cur_idx[drop]] = add;
    cur_idx[add] = cur_idx[drop];
    cenFlag[drop] = 0;
    cenFlag[add]  = 1;

    for (int i = 0; i < cen_len[drop]; i++) {
        int nd = node_of_cen[drop][i];
        node_of_cen[add][i]  = nd;
        idx_of_cen[add][nd]  = i;
        H[nd][add]  = 1;
        H[nd][drop] = 0;
        int pos = idx_of_node[nd][drop];
        cen_of_node[nd][pos]  = add;
        idx_of_node[nd][add]  = pos;
    }
    cen_len[add]  = cen_len[drop];
    cen_len[drop] = 0;
}

void INSERT(int nd, int org_cen, int tar_cen) {
    H[nd][tar_cen] = 1;
    H[nd][org_cen] = 0;

    int pos = idx_of_node[nd][org_cen];
    cen_of_node[nd][pos]     = tar_cen;
    idx_of_node[nd][tar_cen] = pos;

    node_of_cen[tar_cen][cen_len[tar_cen]] = nd;
    idx_of_cen[tar_cen][nd] = cen_len[tar_cen];
    cen_len[tar_cen]++;

    int pos2 = idx_of_cen[org_cen][nd];
    int last = node_of_cen[org_cen][cen_len[org_cen] - 1];
    node_of_cen[org_cen][pos2] = last;
    idx_of_cen[org_cen][last]  = pos2;
    cen_len[org_cen]--;
}

double FLoyd() {
    for (int i = 0; i < verNum; i++)
        for (int j = i; j < verNum; j++)
            compute_edge[i][j] = compute_edge[j][i] = MAX_VAL;

    for (int k = 0; k < p; k++) {
        int tc = cur_cen[k];
        for (int j = 0; j < verNum; j++) {
            if (!cenFlag[j]) {
                for (int jr = 0; jr < r; jr++) {
                    int tc2 = cen_of_node[j][jr];
                    double d = alfa * edge[tc][tc2] + deta * edge[tc2][j];
                    if (d < compute_edge[tc][j]) {
                        compute_edge[tc][j] = d;
                        compute_edge[j][tc] = d;
                    }
                }
            } else {
                compute_edge[tc][j] = compute_edge[j][tc] = alfa * edge[tc][j];
            }
        }
    }

    for (int i = 0; i < verNum; i++) {
        if (cenFlag[i]) continue;
        for (int j = 0; j < verNum; j++) {
            for (int k = 0; k < r; k++) {
                int tc = cen_of_node[i][k];
                double d = gama * edge[i][tc] + compute_edge[tc][j];
                if (d < compute_edge[i][j])
                    compute_edge[i][j] = d;
            }
        }
    }

    double maxVal = 0;
    for (int i = 0; i < verNum; i++)
        for (int j = 0; j < verNum; j++)
            if (compute_edge[i][j] > maxVal)
                maxVal = compute_edge[i][j];

    return maxVal;
}

void init_construct() {
    int clustNum = 0;
    while (clustNum < p) {
        int seed = random_int(verNum);
        while (cenFlag[seed]) seed = random_int(verNum);
        cur_cen[clustNum] = seed;
        cur_idx[seed]     = clustNum;
        cenFlag[seed]     = 1;
        clustNum++;
    }

    for (int i = 0; i < verNum; i++) {
        int u = random_int(p);
        int rrr = 0;
        for (int j = 0; j < p; j++) {
            int tc = cur_cen[(j + u) % p];
            if (tc == i) continue;
            node_of_cen[tc][cen_len[tc]] = i;
            idx_of_cen[tc][i]  = cen_len[tc];
            cen_len[tc]++;
            cen_of_node[i][rrr] = tc;
            idx_of_node[i][tc]  = rrr;
            H[i][tc] = 1;
            rrr++;
            if (rrr == r) break;
        }
    }
}

bool na_first_improve() {
    for (int nd = 0; nd < verNum; nd++) {
        if (cenFlag[nd]) continue;
        for (int ri = 0; ri < r; ri++) {
            int org_cen = cen_of_node[nd][ri];
            for (int k = 0; k < p; k++) {
                int tar_cen = cur_cen[k];
                if (tar_cen == org_cen) continue;
                if (H[nd][tar_cen]) continue;

                INSERT(nd, org_cen, tar_cen);
                double new_val = FLoyd();
                if (new_val < Sc - EPCION) {
                    Sc = new_val;
                    return true;
                } else {
                    INSERT(nd, tar_cen, org_cen);
                }
            }
        }
    }
    return false;
}

bool nh_first_improve() {
    for (int i = 0; i < verNum; i++) U[i] = i;
    int x = verNum;

    while (x > 0) {
        int r0 = random_int(x);
        int add_hub = U[r0];
        x--;
        U[r0] = U[x];

        if (cenFlag[add_hub]) continue;

        for (int j = 0; j < p; j++) {
            int drop_hub = cur_cen[j];
            SWAP(add_hub, drop_hub);
            double new_val = FLoyd();
            if (new_val < Sc - EPCION) {
                Sc = new_val;
                if (Sc < Sb - EPCION) Sb = Sc;
                return true;
            } else {
                SWAP(drop_hub, add_hub);
            }
        }
    }
    return false;
}

void seq_vnd() {
    bool improved = true;
    while (improved) {
        improved = false;
        if (nh_first_improve()) { improved = true; continue; }
        if (na_first_improve())   improved = true;
    }
}

void shake(int k) {
    for (int i = 0; i < k; i++) {
        int add_hub = random_int(verNum);
        while (cenFlag[add_hub]) add_hub = random_int(verNum);
        int drop_hub = cur_cen[random_int(p)];
        SWAP(add_hub, drop_hub);
    }
    Sc = FLoyd();
}

void save_hubs(int* buf) {
    for (int i = 0; i < p; i++) buf[i] = cur_cen[i];
}

void restore_hubs(int* buf) {
    int diff1[100], diff2[100];
    int d1 = 0, d2 = 0;

    for (int i = 0; i < p; i++) {
        bool found = false;
        for (int j = 0; j < p; j++)
            if (buf[i] == cur_cen[j]) { found = true; break; }
        if (!found) diff1[d1++] = buf[i];
    }

    for (int i = 0; i < p; i++) {
        bool found = false;
        for (int j = 0; j < p; j++)
            if (cur_cen[i] == buf[j]) { found = true; break; }
        if (!found) diff2[d2++] = cur_cen[i];
    }

    for (int i = 0; i < d1; i++) SWAP(diff1[i], diff2[i]);
    Sc = FLoyd();
}

void gvns_seq() {
    init_construct();
    Sc = FLoyd();
    Sb = Sc;

    int kmax = p;
    int saved_hubs[100];

    while ((double)(clock() - startTime) / CLOCKS_PER_SEC < cuttingTime) {
        seq_vnd();

        if (Sc < Sb - EPCION) Sb = Sc;

        save_hubs(saved_hubs);
        double saved_Sc = Sc;

        int k = 1;
        while (k <= kmax && (double)(clock() - startTime) / CLOCKS_PER_SEC < cuttingTime) {
            restore_hubs(saved_hubs);
            Sc = saved_Sc;

            shake(k);
            seq_vnd();

            if (Sc < saved_Sc - EPCION) {
                save_hubs(saved_hubs);
                saved_Sc = Sc;
                if (Sc < Sb - EPCION) Sb = Sc;
                k = 1;
            } else {
                k++;
            }
        }

        restore_hubs(saved_hubs);
        Sc = saved_Sc;
    }
}

void write_outcome() {
    ostringstream outfilename;
    outfilename << "p" << p << "r" << r << ".txt";
    ofstream outFile(outfilename.str(), ios_base::app);
    if (!outFile) {
        cerr << "无法打开文件" << endl;
        exit(0);
    }
    outFile << "p:" << p << " r:" << r << " Sb:" << fixed << setprecision(2) << Sb << endl;
    outFile.close();
}

void freememo() {
    for (int i = 0; i < verNum; i++) {
        delete[] edge[i];
        delete[] compute_edge[i];
        delete[] H[i];
        delete[] cen_of_node[i];
        delete[] idx_of_node[i];
        delete[] node_of_cen[i];
        delete[] idx_of_cen[i];
    }
    delete[] edge;
    delete[] compute_edge;
    delete[] H;
    delete[] cur_cen;
    delete[] cenFlag;
    delete[] cur_idx;
    delete[] cen_of_node;
    delete[] idx_of_node;
    delete[] node_of_cen;
    delete[] idx_of_cen;
    delete[] cen_len;
    delete[] U;
}

int main(int argc, char** argv) {
    if (argc > 1) inFile = argv[1];
    if (argc > 2) p  = atoi(argv[2]);
    if (argc > 3) r  = atoi(argv[3]);
    if (argc > 4) cuttingTime = atoi(argv[4]);
    if (argc > 5) gama = atof(argv[5]);
    if (argc > 6) alfa = atof(argv[6]);
    if (argc > 7) deta = atof(argv[7]);

    srand((unsigned)time(NULL));
    read_initial(inFile);

    cout << "File:" << inFile << "  p=" << p << "  r=" << r
        << "  gama=" << gama << "  alfa=" << alfa << "  deta=" << deta << endl;

    startTime = clock();
    gvns_seq();

    cout << fixed << setprecision(2) << "Sb: " << Sb << endl;

    write_outcome();
    freememo();

    return 0;
}
