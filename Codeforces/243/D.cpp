#include <iostream>
#include <cstdio>
#include <string>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <vector>
#include <sstream>
#include <set>
#include <ctime>
#include <map>
#include <stack>
#include <cmath>
using namespace std;
const int N = 100005;
const double pi = acos (-1.0);
const int inf = 1 << 29;
const int lim = 500;
vector <int> a[N];
int n , x[N] , y[N];

bool exist (int x , int y) {
    if (x >= N || y >= N) return false;
    return binary_search (a[x].begin () , a[x].end () , y);
}

int main(){
    cin >> n;
    for (int i = 0 ; i < n ; i ++) {
        cin >> x[i] >> y[i];
        a[x[i]].push_back (y[i]);
    }
    for (int i = 0 ; i < N ; i ++) {
        sort (a[i].begin () , a[i].end ());
    }
    int ans = 0;
    for (int i = 0 ; i < N ; i ++) {
        if (a[i].size() < lim) {
            for (int j = 0 ; j < a[i].size() ; j ++) {
                for (int k = j + 1 ; k < a[i].size() ; k ++) {
                    int d = a[i][k] - a[i][j];
                    if (d > 0 && exist (i + d , a[i][k]) && exist (i + d , a[i][j]))
                        ans ++;
                }
            }
        }
        else {
            for (int j = i + 1 ; j < N ; j ++) {
                for (int k = 0 ; k < a[j].size() ; k ++) {
                    int d = j - i;
                    if (exist (i , a[j][k]) && exist (i , a[j][k] + d) && exist (j , a[j][k] + d))
                        ans ++;
                }
            }
        }
    }
    cout << ans << endl;
    return 0;
}
