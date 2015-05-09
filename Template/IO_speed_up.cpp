//一整行数的输入
char buf[MAXN];
gets(buf);
int v;
char *p=strtok(buf," ");
while(p) {
    sscanf(p,"%d",&v);
    p=strtok(NULL," ");
}

//ACMonster's IO
int get() {
    char c;
    while(c=getchar(),(c<'0'||c>'9')&&(c!='-'));
    bool flag=(c=='-');
    if(flag)
        c=getchar();
    int x=0;
    while(c>='0'&&c<='9') {
        x=x*10+c-'0';
        c=getchar();
    }
    return flag?-x:x;
}
void output(long long x) {
    if(x<0) {
        putchar('-');
        x=-x;
    }
    int len=0,data[20];
    while(x) {
        data[len++]=x%10;
        x/=10;
    }
    if(!len)
        data[len++]=0;
    while(len--)
        putchar(data[len]+'0');
    putchar('\n');
}

inline int readint() {
    char c=getchar();
    while(!isdigit(c))
        c=getchar();
    int x=0;
    while(isdigit(c)) {
        x=x*10+c-'0';
        c=getchar();
    }
    return x;
}
char buf[20];
inline void writeint(int x) {
    if(x==0) {
        putchar('0');
        return;
    }
    if(x<0) {
        putchar('-');
        x=-x;
    }
    int bas=0;
    while(x) {
        buf[bas++]=x%10+'0';
        x/=10;
    }
    while(k--)
        putchar(buf[bas]);
}

template<class T> inline T& RDD(T &x) {
    char c;
    for (c = getchar(); c < '-'; c = getchar());
    if (c == '-') {
        x = '0' - getchar();
        for (c = getchar(); '0' <= c && c <= '9'; c = getchar()) x = x * 10 + '0' - c;
    } else {
        x = c - '0';
        for (c = getchar(); '0' <= c && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    }
    return x;
}
