#include <cstdio>
#include <algorithm>
using namespace std;

typedef long long LL;
#define lson l , m , rt << 1
#define rson m + 1 , r , rt << 1 | 1
const int maxn = 111111;
LL sum[maxn<<2];
int MAX[maxn<<2];
void PushUP(int rt) {
    MAX[rt] = max(MAX[rt<<1],MAX[rt<<1|1]);
	sum[rt] = sum[rt<<1] + sum[rt<<1|1];
}
void build(int l,int r,int rt) {
	if (l == r) {
		scanf("%I64d",&sum[rt]);
        MAX[rt] = sum[rt];
		return ;
	}
	int m = (l + r) >> 1;
	build(lson);
	build(rson);
	PushUP(rt);
}
void mody(int L,int R,int mod,int l,int r,int rt){
    if(MAX[rt] < mod)return;
    if(l == r){
        MAX[rt] = sum[rt] = sum[rt] % mod;
        return;
    }
    int m = (l + r) >> 1;
    if(L <= m)mody(L, R, mod, lson);
    if(R > m) mody(L, R, mod, rson);
    PushUP(rt);
}

void modify(int p,int val,int l,int r,int rt) {
	if (l == r) {
		MAX[rt] = sum[rt] = val;
		return ;
	}
	int m = (l + r) >> 1;
	if (p <= m) modify(p , val , lson);
	else modify(p , val , rson);
	PushUP(rt);
}
LL query(int L,int R,int l,int r,int rt) {
	if (L <= l && r <= R) {
		return sum[rt];
	}
	int m = (l + r) >> 1;
	LL ret = 0;
	if (L <= m) ret += query(L , R , lson);
	if (R > m) ret += query(L , R , rson);
	return ret;
}
int n, m;
int main()
{
    scanf("%d%d",&n,&m);
    build(1, n, 1);
    for(int i = 0 ;i < m ; i++){
        int op;
        scanf("%d",&op);
        if(op == 1){
            int l,r;
            scanf("%d%d",&l,&r);
            printf("%I64d\n",query(l, r, 1, n, 1));
        }else if(op == 2){
            int l,r,x;
            scanf("%d%d%d",&l,&r,&x);
            mody(l, r, x, 1, n, 1);
        }else if(op == 3){
            int k,x;
            scanf("%d%d",&k,&x);
            modify(k, x, 1, n, 1);
        }
    }
    return 0;
}


