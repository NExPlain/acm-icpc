#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<vector>
#include<algorithm>
using namespace std;

const int N = 200000+20;
typedef pair<int,int> PII;
#define CLR(a,b) memset(a,b,sizeof(a))
bool check[N];
int a[N];
int pos[N];
int n;
int prime[N];
int tot;

void shai(int N)
{
    CLR(check,0);
    tot = 0;
    for(int i = 2; i <= N; i ++){
        if(!check[i]){
            prime[tot++] = i;
        }
        for(int j = 0 ; j < tot ; j++){
            if(i * prime[j] > N)break;
            check[i * prime[j]] = true;
            if(i % prime[j] == 0)break;
        }
    }
}

void solve()
{
    vector<PII> ans;
    for(int i = 1; i <= n ; i++){
        while(pos[i] != i){
            int cha = pos[i] - i + 1;
            int wei = upper_bound(prime, prime+tot, cha) - prime;
            wei--;
            ans.push_back(make_pair(pos[i]-prime[wei]+1,pos[i]));
            int t1 = pos[i]-prime[wei]+1,t2 = pos[i];
            swap(a[t1],a[t2]);
            swap(pos[a[t1]],pos[a[t2]]);
        }
    }
    cout << ans.size() << endl;
    for(int i = 0 ;i < ans.size(); i ++){
        cout << ans[i].first << " " << ans[i].second << endl;
    }
}
int main()
{
    shai(N-1);
    cin >> n ;
    for(int i = 1 ; i <= n ; i++){
        cin >> a[i];
        pos[a[i]] = i;
    }
    solve();
    return 0;
}
