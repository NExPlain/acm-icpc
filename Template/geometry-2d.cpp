#include<cstdio>
#include<cmath>
#include<algorithm>
#include<vector>
using namespace std;
const int maxn=100;
const double INF=1e100;
const double eps=1e-8;
const double pi=acos(-1.0);
#define next(i) ((i+1)%n)
inline int cmp(double x) { //判断数的符号
    return x<-eps?-1:x>eps?1:0;
}
inline double sqr(double x) { //计算数的平方
    return x*x;
}
struct point { //点类
    double x,y;
    point() {}
    point(double a,double b):x(a),y(b) {}
    void input() {
        scanf("%lf%lf",&x,&y);
    }
    friend point operator+(const point &a,const point &b) {
        return point(a.x+b.x,a.y+b.y);
    }
    friend point operator-(const point &a,const point &b) {
        return point(a.x-b.x,a.y-b.y);
    }
    friend bool operator==(const point &a,const point &b) {
        return cmp(a.x-b.x)==0&&cmp(a.y-b.y)==0;
    }
    friend point operator*(const point &a,const double b) {
        return point(a.x*b,a.y*b);
    }
    friend point operator*(const double a,const point &b) {
        return point(a*b.x,a*b.y);
    }
    friend point operator/(const point &a,const double b) {
        return point(a.x/b,a.y/b);
    }
    double norm() {
        return sqrt(sqr(x)+sqr(y));
    }
};
double det(const point &a,const point &b) { //计算向量叉积
    return a.x*b.y-a.y*b.x;
}
double dot(const point &a,const point &b) { //计算向量点积
    return a.x*b.x+a.y*b.y;
}
double dist(const point &a,const point &b) { //计算两点距离
    return (a-b).norm();
}
point rotate_point(const point &p,double A) { //计算向量p绕逆时针旋转A弧度得到的向量
    return point(p.x*cos(A)-p.y*sin(A),p.x*sin(A)+p.y*cos(A));
}
struct line { //线段或直线类
    point a,b;
    line() {}
    line(point x,point y):a(x),b(y) {}
};
line point_make_line(const point &a,const point &b) { //生成线段或直线
    return line(a,b);
}
double dis_point_segment(const point &p,const point &s,const point &t) { //计算点p到线段st距离
    if(cmp(dot(p-s,t-s))<0)
        return (p-s).norm();
    if(cmp(dot(p-t,s-t))<0)
        return (p-t).norm();
    return fabs(det(s-p,t-p)/dist(s,t));
}

void PointProjLine(const point &p,const point &s,const point &t,point &cp) { //计算点p到线段st的垂足，保存于cp
    double r=dot(t-s,p-s)/dot(t-s,t-s);
    cp=s+r*(t-s);
}
bool PointOnSegment(const point &p,const point &s,const point &t) { //判断点p是否在线段sp上
    return cmp(det(p-s,t-s))==0&&cmp(dot(p-s,p-t))<=0;
}
bool parallel(const line &a,const line &b) { //判断直线a和b是否平行
    return !cmp(det(a.a-a.b,b.a-b.b));
}
bool line_make_point(const line &a,const line &b,point &ret) { //判断直线a和b是否相交，若相交交点保存于ret
    if(parallel(a,b))
        return false;
    double s1=det(a.a-b.a,b.b-b.a);
    double s2=det(a.b-b.a,b.b-b.a);
    ret=(s1*a.b-s2*a.a)/(s1-s2);
    return true;
}
line move_d(const line &a,const double len) { //计算直线a沿法向量平移距离len得到的直线
    point d=a.b-a.a;
    d=d/d.norm();
    d=rotate_point(d,pi/2);
    return line(a.a+d*len,a.b+d*len);
}
double direction(const point &a,const point &b,const point &c) {
    return det(c-a,b-a);
}
bool OnSegment(const point &a,const point &b,const point &c) {
    if(c.x>=min(a.x,b.x)&&c.x<=max(a.x,b.x)&&c.y>=min(a.y,b.y)&&c.y<=max(a.y,b.y))
        return true;
    return false;
}
bool SegmentsIntersect(const line &a,const line &b) { //判断线段a和b是否相交
    double d1=direction(b.a,b.b,a.a);
    double d2=direction(b.a,b.b,a.b);
    double d3=direction(a.a,a.b,b.a);
    double d4=direction(a.a,a.b,b.b);
    if(d1*d2<0&&d3*d4<0)
        return true;
    else if(d1==0&&OnSegment(b.a,b.b,a.a))
        return true;
    else if(d2==0&&OnSegment(b.a,b.b,a.b))
        return true;
    else if(d3==0&&OnSegment(a.a,a.b,b.a))
        return true;
    else if(d4==0&&OnSegment(a.a,a.b,b.b))
        return true;
    return false;
}
struct polygon { //多边形类
    int n;
    point a[maxn];
    polygon() {}
    double perimeter() { //计算多边形周长
        double sum=0;
        a[n]=a[0];
        for(int i=0; i<n; ++i)
            sum+=(a[i+1]-a[i]).norm();
        return sum;
    }
    double area() { //计算多边形面积
        double sum=0;
        a[n]=a[0];
        for(int i=0; i<n; ++i)
            sum+=det(a[i+1],a[i]);
        return sum/2;
    }
    int Point_In(const point &t) { //判断点t是否在多边形内部，0表示在外，1表示在内，2表示在边界,O(n)
        int num=0,d1,d2,k;
        a[n]=a[0];
        for(int i=0; i<n; ++i) {
            if(PointOnSegment(t,a[i],a[i+1]))
                return 2;
            k=cmp(det(a[i+1]-a[i],t-a[i]));
            d1=cmp(a[i].y-t.y);
            d2=cmp(a[i+1].y-t.y);
            if(k>0&&d1<=0&&d2>0)
                ++num;
            if(k<0&&d2<=0&&d1>0)
                --num;
        }
        return num!=0;
    }
    point MassCenter() { //计算多边形重心,O(n)
        point ans=point(0,0);
        if(cmp(area())==0) //多边形面积为0时重心无定义，需特殊处理
            return ans;
        a[n]=a[0];
        for(int i=0; i<n; ++i)
            ans=ans+(a[i]+a[i+1])*det(a[i+1],a[i]);
        return ans/area()/6;
    }
    //以下函数需要多边形顶点为整点
    int Border_Int_Point_Num() { //计算多边形边界格点（整点）个数
        int num=0;
        a[n]=a[0];
        for(int i=0; i<n; ++i)
            num+=__gcd(abs(int(a[i+1].x-a[i].x)),abs(int(a[i+1].y-a[i].y)));
        return num;
    }
    int Inside_Int_Point_Num() { //计算多边形内部格点（整点）个数
        return int(area())+1-Border_Int_Point_Num()/2;
    }
};
struct polygon_convex { //凸多边形类
    vector<point> P;
    polygon_convex(int Size=0) {
        P.resize(Size);
    }
};
bool comp_less(const point &a,const point &b) {
    return cmp(a.x-b.x)<0||cmp(a.x-b.x)==0&&cmp(a.y-b.y)<0;
}
polygon_convex convex_hull(vector<point> &a) { //计算点集a的凸包，逆时针序排列，水平序求法以避免精度问题，O(nlogn)
    polygon_convex ret(2*a.size()+5);
    sort(a.begin(),a.end(),comp_less);
    a.erase(unique(a.begin(),a.end()),a.end());
    int m=0;
    for(int i=0; i<a.size(); ++i) {
        while(m>1&&cmp(det(ret.P[m-1]-ret.P[m-2],a[i]-ret.P[m-2]))<=0)
            --m;
        ret.P[m++]=a[i];
    }
    int k=m;
    for(int i=int(a.size())-2; i>=0; --i) {
        while(m>k&&cmp(det(ret.P[m-1]-ret.P[m-2],a[i]-ret.P[m-2]))<=0)
            --m;
        ret.P[m++]=a[i];
    }
    ret.P.resize(m);
    if(a.size()>1)
        ret.P.resize(m-1);
    return ret;
}
bool containOn(const polygon_convex &a,const point &b) { //判断点b是否在凸包a内部，true表示在内或边界，O(n)
    int n=a.P.size();
    int sign=0;
    for(int i=0; i<n; ++i) {
        int x=cmp(det(a.P[i]-b,a.P[next(i)]-b));
        if(x) {
            if(sign) {
                if(sign!=x)
                    return false;
            } else
                sign=x;
        }
    }
    return true;
}
int containOlogn(const polygon_convex &a,const point &b) { //判断点b是否在凸包a内部，true表示在内或边界，O(logn)
    int n=a.P.size();
    point g=(a.P[0]+a.P[n/3]+a.P[2*n/3])/3; //找一个凸包内部的点g
    int l=0,r=n;
    while(l+1<r) { //二分凸包 g-a.P[a]-a.P[b]
        int mid=(l+r)/2;
        if(cmp(det(a.P[l]-g,a.P[mid]-g))>0) {
            if(cmp(det(a.P[l]-g,b-g))>=0&&cmp(det(a.P[mid]-g,b-g))<0)
                r=mid;
            else
                l=mid;
        } else {
            if(cmp(det(a.P[l]-g,b-g))<0&&cmp(det(a.P[mid]-g,b-g))>=0)
                l=mid;
            else
                r=mid;
        }
    }
    r%=n;
    int z=cmp(det(a.P[r]-b,a.P[l]-b))-1;
    if(z==-2)
        return 1;
    return z;
}
double convex_diameter(polygon_convex &a,int &First,int &Second) { //计算凸包a上最远欧几里得距离，最远点对应标号保存于First和Second，O(n)
    vector<point> &p=a.P;
    int n=p.size();
    double maxd=0;
    if(n==1) {
        First=Second=0;
        return maxd;
    }
    for(int i=0,j=1; i<n; ++i) {
        while(cmp(det(p[next(i)]-p[i],p[j]-p[i])-det(p[next(i)]-p[i],p[next(j)]-p[i]))<0)
            j=next(j);
        double d=dist(p[i],p[j]);
        if(d>maxd) {
            maxd=d;
            First=i;
            Second=j;
        }
        d=dist(p[next(i)],p[next(j)]);
        if(d>maxd) {
            maxd=d;
            First=i;
            Second=j;
        }
    }
    return maxd;
}
double convex_perimeter(polygon_convex &a) { //计算凸多边形周长
    double sum=0;
    for(int i=1; i<a.P.size(); ++i)
        sum+=(a.P[i]-a.P[i-1]).norm();
    sum+=(a.P[0]-a.P[a.P.size()-1]).norm();
    return sum;
}
double convex_area(polygon_convex &a) { //计算凸多边形面积
    double sum=0;
    for(int i=1; i<a.P.size(); ++i)
        sum+=det(a.P[i],a.P[i-1]);
    sum+=det(a.P[0],a.P[a.P.size()-1]);
    return sum/2;
}
struct halfPlane {
    //ax+by+c<=0;即向量左侧平面
    double a,b,c;
    halfPlane(const point &p,const point &q) {
        a=p.y-q.y;
        b=q.x-p.x;
        c=det(p,q);
    }
    halfPlane(double aa,double bb,double cc):a(aa),b(bb),c(cc) {}
};
//计算点a带入到直线方程中的函数值
double calc(const halfPlane &L,const point &a) {
    return a.x*L.a+a.y*L.b+L.c;
}
//求点a和b连线与半平面L的交点
point Intersect(const point &a,const point &b,const halfPlane &L) {
    point ret;
    double t1=calc(L,a),t2=calc(L,b);
    ret.x=(t2*a.x-t1*b.x)/(t2-t1);
    ret.y=(t2*a.y-t1*b.y)/(t2-t1);
    return ret;
}
//将一个凸多边形和一个半平面交
polygon_convex cut(const polygon_convex &a,const halfPlane &L) {
    int n=a.P.size();
    polygon_convex ret;
    for(int i=0; i<n; ++i)
        if(calc(L,a.P[i])<-eps)
            ret.P.push_back(a.P[i]);
        else {
            int j=i-1;
            if(j<0)
                j=n-1;
            if(calc(L,a.P[i])<-eps)
                ret.P.push_back(Intersect(a.P[j],a.P[i],L));
            j=i+1;
            if(j==n)
                j=0;
            if(calc(L,a.P[j])<-eps)
                ret.P.push_back(Intersect(a.P[i],a.P[j],L));
        }
    return ret;
}
//求一个多边形的核
polygon_convex core(polygon &a) {
    polygon_convex ret;
    ret.P.push_back(point(-INF,-INF));
    ret.P.push_back(point(INF,-INF));
    ret.P.push_back(point(INF,INF));
    ret.P.push_back(point(-INF,INF));
    int n=a.n;
    for(int i=0; i<n; ++i) {
        halfPlane L(a.a[i],a.a[(i+1)%n]);
        ret=cut(ret,L);
    }
    return ret;
}
double mysqrt(double n) {
    return sqrt(max(0.0,n));
}
//求圆与线段（直线）的交点
void circle_cross_line(const point &a,const point &b,const point &o,double r,point ret[],int &num) {
    double x0=o.x,y0=o.y,x1=a.x,y1=a.y,x2=b.x,y2=b.y;
    double dx=x2-x1,dy=y2-y1;
    double A=sqr(dx)+sqr(dy),B=2*dx*(x1-x0)+2*dy*(y1-y0),C=sqr(x1-x0)+sqr(y1-y0)-sqr(r);
    double delta=sqr(B)-4*A*C;
    num=0;
    if(cmp(delta)>=0) {
        double t1=(-B-mysqrt(delta))/(2*A);
        double t2=(-B+mysqrt(delta))/(2*A);
        //去掉以下两个if判断，即可判断与直线交点
        if(cmp(t1-1)<=0&&cmp(t1)>=0)
            ret[num++]=point(x1+t1*dx,y1+t1*dy);
        if(cmp(t2-1)<=0&&cmp(t2)>=0)
            ret[num++]=point(x1+t2*dx,y1+t2*dy);
    }
}
//求圆与简单多边形交的面积，圆处于原点
point crosspt(const point &a,const point &b,const point &p,const point &q) {
    double a1=det(b-a,p-a);
    double a2=det(b-a,q-a);
    return (p*a2-q*a1)/(a1-a2);
}
double sector_area(const point &a,const point &b,double r) {
    double theta=atan2(a.y,a.x)-atan2(b.y,b.x);
    while(theta<=0)
        theta+=2*pi;
    while(theta>2*pi)
        theta-=2*pi;
    theta=min(theta,2*pi-theta);
    return sqr(r)*theta/2;
};
double calc(point &a,point &b,double r) {
    point p[2];
    int num=0,ina=cmp(a.norm()-r)<0,inb=cmp(b.norm()-r)<0;
    if(ina) {
        if(inb)
            return fabs(det(a,b))/2.0;
        else {
            circle_cross_line(a,b,point(0,0),r,p,num);
            return sector_area(b,p[0],r)+fabs(det(a,p[0]))/2.0;
        }
    } else {
        circle_cross_line(a,b,point(0,0),r,p,num);
        if(inb)
            return sector_area(p[0],a,r)+fabs(det(p[0],b))/2.0;
        else {
            if(num==2)
                return sector_area(a,p[0],r)+sector_area(p[1],b,r)+fabs(det(p[0],p[1]))/2.0;
            else
                return sector_area(a,b,r);
        }
    }
}
double area(point res[],double r,int n) { //需保证res[n]==res[0]
    double ret=0;
    for(int i=0; i<n; ++i) {
        int sgn=cmp(det(res[i],res[i+1]));
        if(sgn)
            ret+=sgn*calc(res[i],res[i+1],r);
    }
    return ret;
}
//求最小圆覆盖（需要先把a[]中的点打乱）
void circle_center(const point &p0,const point &p1,const point &p2,point cp) {
    double a1=p1.x-p0.x,b1=p1.y-p0.y,c1=(sqr(a1)+sqr(b1))/2;
    double a2=p2.x-p0.x,b2=p2.y-p0.y,c2=(sqr(a2)+sqr(b2))/2;
    double d=a1*b2-a2*b1;
    cp.x=p0.x+(c1*b2-c2*b1)/d;
    cp.y=p0.y+(a1*c2-a2*c1)/d;
}
void circle_center(const point &p0,const point &p1,point &cp) {
    cp.x=(p0.x+p1.x)/2;
    cp.y=(p0.y+p1.y)/2;
}
bool point_in(const point &p,const point &o,double r) {
    return cmp(dist(p,o)-r)<0;
}
void min_circle_cover(const point a[],int n,point &o,double &r) {
    r=0;
    o=a[0];
    for(int i=1; i<n; ++i)
        if(!point_in(a[i],o,r)) {
            o=a[i];
            r=0;
            for(int j=0; j<i; ++j)
                if(!point_in(a[j],o,r)) {
                    circle_center(a[i],a[j],o);
                    r=dist(a[j],o);
                    for(int k=0; k<j; ++k)
                        if(!point_in(a[k],o,r)) {
                            circle_center(a[i],a[j],a[k],o);
                            r=dist(a[k],o);
                        }
                }
        }
}
//求两个圆的交点
point rotate(const point &p,double cost,double sint) {
    double x=p.x,y=p.y;
    return point(x*cost-y*sint,x*sint+y*cost);
}
pair<point,point> crosspoint(const point &ap,double ar,const point &bp,double br) {
    double d=dist(ap,bp),cost=(sqr(ar)+sqr(d)-sqr(br))/(2*ar*d);
    double sint=sqrt(1.0-sqr(cost));
    point v=(bp-ap)/dist(bp-ap)*ar;
    return make_pair(ap+rotate(v,cost,-sint),ap+rotate(v,cost,sint));
}
int main() {
    return 0;
}
