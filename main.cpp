#include <iostream>
#include <vector>
#include <random>
#include <queue>
#include <utility>
#include <functional> // greater/less (comparison for priority queue)
#include <time.h>
#include <algorithm>

using namespace std;


class Obj
{
public:
    Obj(int pX, int pY) : x(pX), y(pY) { cluster = 0; }
    double dist(Obj& to) { return (x - to.getX()) * (x - to.getX()) + (y - to.getY())*(y - to.getY()); }
    int getX() { return x; }
    int getY() { return y; }
    int getCluster() { return cluster; }

private:
    int x;
    int y;
    int cluster;
};

double distSq(int x1, int x2, int y1, int y2) {
    return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
}
void printPair(pair<double, int> par) {
    cout << par.first << ", " << par.second << endl;
}

void flame(vector<Obj> objects)
{
    for( Obj o : objects ) {
        cout << o.getX() << " " << o.getY() << endl;
    }
}

void flame(vector<int> p_xs, vector<int> p_ys)
{
    for( int i = 0; i < p_xs.size(); i++ ) {
        cout << p_xs[i] << " " << p_ys[i] << endl;
    }
}





// k nearest neighbours
int* flame(int* p_xs, int* p_ys, int size, int k = 0)
{
    if ( k == 0 )
        k = min(40, size / 5);
    cout << "flame started\n";

    int** indx = new int*[size];
    double** dist = new double*[size];
    double* density = new double[size];

    auto start = time(0);

    // should be O(n^2 log(k)) time ( n times repeat: insert n elements into container of size k with insertion log(k) )
    for( int i = 0; i < size; i++ ) {
        indx[i] = new int[k];
        dist[i] = new double[k];
//        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double,int>> > myQ;
        priority_queue< pair<double, int> > myQ;
        for( int j = 0; j < size; j++ ) {
            if( j == i ) continue;
            myQ.push(make_pair(distSq(p_xs[i], p_xs[j], p_ys[i], p_ys[j]) , j));

            // largest are first, so we pop the right one
            if( myQ.size() > k )
                myQ.pop();
        }

        density[i] = 0;
        for(int x = 0; x < k; x++) {
            auto par = myQ.top();
            myQ.pop();
            dist[i][x] = par.first;
            indx[i][x] = par.second;
            density[i] += par.first;

//            printPair(par);
        }
    }
    auto end = time(0);
    cout << "time taken: " << difftime(end, start) << endl << endl;
    cout << "distances and k-nn calculated" << endl;


    // ok maybe not best and worst, but unchangable, and setup membership here already !!! aja ne, ka ne vemo kolko clusterjev je !!!

    // probi se z enim prehodom, tko da mas counter, kolko jih je vecji, ce 0 -> neki, ce vsi -> neki druzga
    bool* best = new bool[size];
    memset(best, true, size);
    int numClusters = size;
    for( int i = 0; i < size; i++ ) {
        for( int nn = 0; nn < k; nn++ ) {
            if( density[indx[i][nn]] > density[i] ) {
                best[i] = false;
                numClusters--;
                break;
            }
            // cant do this, cause if x in knn of y =/=> y in knn x
//            else best[nn] = false;
        }
    }
    cout << "highest density done" << endl;
    cout << "num clusters " << numClusters << endl;

    bool* worst = new bool[size];
    memset(worst, true, size);
    for( int i = 0; i < size; i++ ) {
        for( int nn = 0; nn < k; nn++ ) {
            if( density[indx[i][nn]] < density[i] ) {
                worst[i] = false;
                break;
            }
        }
    }
    cout << "lowest density done" << endl;

    // outliers cluster
    numClusters++;


    double** membership = new double*[size];
    int bestIndx = 0;
    for( int i = 0; i < size; i++ )
    {
        membership[i] = new double[numClusters];
        if( worst[i] ) {
            memset(membership[i], 0, numClusters  * sizeof(double));
            membership[i][numClusters-1] = 1;
            best[i] = true; // ok, cause else if... we do this, so we have a vector of weather a membership is fixed or not
        }
        else if( best[i] ) {
//            cout << "best index" << bestIndx << endl;
            fill( membership[i], membership[i] + numClusters , 0 );
            membership[i][bestIndx++] = 1;
        }
        else {
            fill( membership[i], membership[i] + numClusters , 1.0 / numClusters);
        }
    }
    cout << "membership initialised\n";


    // we change the dists to the weights
    double totalDist = 0;
    double totalDistN = 0;
    for( int i = 0; i < size; i++ ) {
        totalDist = 0;
        for( int j = 0; j < k; j++ ) {
            dist[i][j] = sqrt(dist[i][j]);
            totalDist += dist[i][j];
        }
        totalDistN = totalDist * (k-1);
        cout << dist[i][0] << " " << totalDist << " ";
        for( int j = 0; j < k; j++ ) {
            dist[i][j] = (totalDist - dist[i][j]) / totalDistN;
//            cout << dist[i][j] << ",";
        }
        cout << dist[i][0] << endl;
//        cout << endl;
    }

    double** membership2 = new double*[size];

    for( int i = 0; i < size; i++ ) {
        membership2[i] = new double[numClusters];
    }

    for(int iter = 0; iter < 100; iter++ ) {
        for( int i = 0; i < size; i++ ) {
            if( best[i] ) continue;

            memset(membership2[i], 0, numClusters * sizeof(double)); // reset to 0
            //fill in new values
            for( int j = 0; j < k; j++ ) {
                for( int m = 0; m < numClusters; m++ ) {
                    membership2[i][m] += (membership[indx[i][j]][m] * dist[i][j] );
                }
            }


//            double total = 0;
//            for( int m = 0; m < numClusters; m++ ) {
//                total += membership[i][m];
//            }
//            cout << total << endl;
        }
        for( int i = 0; i < size; i++ ) {
            if( best[i] ) continue;
//            delete[] membership[i];
            membership[i] = membership2[i];
        }

    }
    cout << "calculated membership\n";

//    for(int i= 0; i < size; i++) {
//        for( int m = 0; m < numClusters; m++ ) {
//            cout << membership[i][m] <<", ";
//        }
//        cout << endl;
//    }

    int* member = new int[size];
    // max membership
    for( int i = 0; i < size; i++ ) {
        int max = 0;
        double maxVal = 0;
        for( int j = 0; j < numClusters; j++ ) {
            if( membership[i][j] > maxVal ) {
                max = j;
                maxVal = membership[i][j];
            }
//            cout << i << " " << membership[i][j]<< endl;
        }
        member[i] = max;
    }

    return member;
}





void test( int n = 100 )
{
    cout << "start test\n";
    int x, y;
    vector<int> v_xs, v_ys;
    vector<Obj> objects;
    int* xs = new int[n];
    int* ys = new int[n];

    // not checking if repeated point
    for( int i = 0; i < n; i++ ) {
        x = rand() % 1000;
        y = rand() % 1000;

        objects.push_back( Obj(x, y) );
        v_xs.push_back(x);
        v_ys.push_back(y);
        xs[i] = x;
        ys[i] = y;
    }

    int * result = flame(xs, ys, n);
    cout << "flame complete" << endl;
    for( int i = 0; i < n; i++ ) {
        cout << result[i] << ", ";
    }
    cout << endl;

    delete[] xs;
    delete[] ys;
    cout << "finished" << endl;
}

int main()
{
    cout << "start main" << endl;
//    Obj a = Obj(21, 2);

//    vector<Obj> objects;
//    objects.push_back(a);

//    flame(objects);

    test(100);

    return 0;
}

