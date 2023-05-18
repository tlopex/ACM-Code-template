const double eps=1e-8;
const double pi=acos(-1.0);
const int maxp=1010;
#define equals(a,b) (fabs((a)-(b))<eps)

const int COUNTER_CLOCKWISE = 1;
const int CLOCKWISE = -1;
const int ONLINE_BACK = 2;
const int ONLINE_FRONT = -2;
const int ON_SEGMENT = 0;

//Intercsect Circle & Circle
const int ICC_SEPERATE = 4;
const int ICC_CIRCUMSCRIBE = 3;
const int ICC_INTERSECT = 2;
const int ICC_INSCRIBE = 1;
const int ICC_CONTAIN = 0;

int sgn(double x)
{
    if(fabs(x)<eps) return 0;
    else return x<0?-1:1;
}
int Dcmp(double x,double y)//比较两个浮点数
{
    if(fabs(x-y)<eps) return 0;
    else return x<y?-1:1;
}
//平面几何  点和线
struct Point
{
    double x,y;
    Point(){};
    Point(double x,double y):x(x),y(y){};
    Point operator + (Point B){return Point(x+B.x,y+B.y);}
    Point operator - (Point B){return Point(x-B.x,y-B.y);}
    Point operator * (double k){return Point(x*k,y*k);}
    Point operator / (double k){return Point(x/k,y/k);}
    bool operator == (Point B){return sgn(x-B.x)==0&&sgn(y-B.y)==0;}
    bool operator < (Point B){return sgn(x-B.x)<0||(sgn(x-B.x)==0&&sgn(y-B.y)<0);}
    double norm(){return x*x+y*y;}
    double abs(){return sqrt(norm());}
};
typedef Point Vector;
typedef vector<Point> Polygon;
double Dot(Vector A,Vector B){return A.x*B.x+A.y*B.y;}//点积
double Len(Vector A){return sqrt(Dot(A,A));}
double Angle(Vector A,Vector B){return acos(Dot(A,B)/Len(A)/Len(B));} //A与B的夹角
double Cross(Vector A,Vector B) {return A.x*B.y-A.y*B.x;}//叉积
double Area2(Point A,Point B,Point C){return Cross(B-A,C-A);}//三角形ABC面积的两倍
double Distance(Point A,Point B){return hypot(A.x-B.x,A.y-B.y);}//hypot是求斜边的 自带函数
double Norm(Vector a){return a.x*a.x+a.y*a.y;}
double Abs(Vector a){return sqrt(Norm(a));}
bool isOrthogonal(Vector a,Vector b){return sgn(Dot(a,b))==0;}
bool isOrthogonal(Point a1,Point a2,Point b1,Point b2){return isOrthogonal(a1-a2,b1-b2);}
bool isOrthogonal(Segment s1,Segment s2){return sgn(Dot(s1.p2-s1.p1,s2.p2-s2.p1))==0;}
bool isParallel(Vector A,Vector B){return sgn(Cross(A,B))==0;}
bool isParallel(Point a1,Point a2,Point b1,Point b2){return isParallel(a1-a2,b1-b2);}
bool isParallel(Segment s1,Segment s2){return sgn(Cross(s1.p2-s1.p1,s2.p2-s2.p1))==0; }
Vector Normal(Vector A){   double L = Abs(A);  return Vector(-A.y / L, A.x / L);}//单位法向量
Vector Rotate(Vector A,double rad)
{
    return Vector(A.x*cos(rad)-A.y*sin(rad),A.x*sin(rad)+A.y*cos(rad));
} //向量A逆时针旋转rad度
struct Line
{
    Point p1,p2;
    Line(){};
    Line(Point p1,Point p2):p1(p1),p2(p2){};
    Line(Point p,double angle)//0<=angle<pi
    {
        p1=p;
        if(sgn(angle-pi/2)==0){p2=(p1+Point(0,1));}
        else{p2=(p1+Point(1,tan(angle)));}
    }
    Line(double a,double b,double c)
    {
        if(sgn(a)==0)
        {
            p1=Point(0,-c/b);
            p2=Point(1,-c/b);
        }
        else if(sgn(b)==0)
        {
            p1=Point(-c/a,0);
            p2=Point(-c/a,1);
        }
        else
        {
            p1=Point(0,-c/b);
            p2=Point(1,(-c-a)/b);
        }
    }
};
typedef Line Segment;

double Line_angle(Line v)
{
    double k=atan2(v.p2.y-v.p1.y,v.p2.x-v.p1.x);
    if(sgn(k)<0) k+=pi;
    if(sgn(k-pi)==0) k-=pi;
    return k;
}//返回直线倾斜角


int Point_line_relation(Point p,Line v)
{
    int c=sgn(Cross(p-v.p1,v.p2-v.p1));
    if(c<0) return 1;
    if(c>0) return 2;
    return 0;
}//点和直线的关系：1在左侧 2在右侧 0在直线上

bool Point_on_seg(Point p, Line v) { // 0为点不在线段v上；1为点在线段v上
    return sgn(Cross(p - v.p1, v.p2 - v.p1)) == 0 && sgn(Dot(p - v.p1, p - v.p2)) <= 0;
}

int Line_relation(Line v1,Line v2)
{
    if(sgn(Cross(v1.p2-v1.p1,v2.p2-v2.p1))==0)
    {
        if(Point_line_relation(v1.p1,v2)==0) return 1;
        else return 0;
    }
    return 2;
}//两直线的关系 0为平行 1为重合 2为相交

double Dis_point_line(Point p,Line v)
{
    return fabs(Cross(p-v.p1,v.p2-v.p1))/Distance(v.p1,v.p2);
}//点到直线的距离

Point project(Point p,Line s){
  Vector base=s.p2-s.p1;
  double r=Dot(p-s.p1,base)/Norm(base);
  return s.p1+base*r;
}//点在直线上的投影

Point symmetry(Point p,Line v)
{
    Point q=project(p,v);
    return Point(2*q.x-p.x,2*q.y-p.y);
} //点p对直线v的对称点

double Dis_point_seg(Point p,Segment v)
{
    if(sgn(Dot(p-v.p1,v.p2-v.p1))<0||sgn(Dot(p-v.p2,v.p1-v.p2))<0)
    return min(Distance(p,v.p1),Distance(p,v.p2));
    //点的投影不在线段上
    return Dis_point_line(p,v);
    //点的投影在线段上
}
Point getCrossPoint(Segment s1,Segment s2){
  Vector base=s2.p2-s2.p1;
  double d1=abs(Cross(base,s1.p1-s2.p1));
  double d2=abs(Cross(base,s1.p2-s2.p1));
  double t=d1/(d1+d2);
  return s1.p1+(s1.p2-s1.p1)*t;
}

Point getCrossPoint(Point p0,Point p1,Point p2,Point p3){
  return getCrossPoint(Segment(p0,p1),Segment(p2,p3));
}

Point getCrossPointLL(Line l1,Line l2){
  double a=Cross(l1.p2-l1.p1,l2.p2-l2.p1);
  double b=Cross(l1.p2-l1.p1,l1.p2-l2.p1);
  if(abs(a)<eps&&abs(b)<eps) return l2.p1;
  return l2.p1+(l2.p2-l2.p1)*(b/a);
}

int ccw(Point p0,Point p1,Point p2){
  Vector a = p1-p0;
  Vector b = p2-p0;
  if(Cross(a,b) > eps) return COUNTER_CLOCKWISE;//1
  if(Cross(a,b) < -eps) return CLOCKWISE;//-1
  if(Dot(a,b) < -eps) return ONLINE_BACK;//2
  if(a.norm()<b.norm()) return ONLINE_FRONT;//-2
  return ON_SEGMENT;//0
}
bool intersect(Point p1,Point p2,Point p3,Point p4){
    //求两线段是否相交
  return (ccw(p1,p2,p3)*ccw(p1,p2,p4) <= 0 &&
      ccw(p3,p4,p1)*ccw(p3,p4,p2) <= 0 );
}

bool intersect(Segment s1,Segment s2){
  return intersect(s1.p1,s1.p2,s2.p1,s2.p2);
}


double getDistanceLP(Line l,Point p){
  return abs(Cross(l.p2-l.p1,p-l.p1)/Abs(l.p2-l.p1));
}

double getDistanceSP(Segment s,Point p){
  if(Dot(s.p2-s.p1,p-s.p1)<0.0) return Abs(p-s.p1);
  if(Dot(s.p1-s.p2,p-s.p2)<0.0) return Abs(p-s.p2);
  return getDistanceLP(s,p);
}

double getDistance(Segment s1,Segment s2){
  if(intersect(s1,s2)) return 0.0;//两线交叉距离为0
  return min(min(getDistanceSP(s1,s2.p1),getDistanceSP(s1,s2.p2)),
         min(getDistanceSP(s2,s1.p1),getDistanceSP(s2,s1.p2)));
}

double getDistance(Point p0,Point p1,Point p2,Point p3){
  return getDistance(Segment(p0,p1),Segment(p2,p3));
}

//凸包
istream &operator >> (istream &is,Point &p){
  is>>p.x>>p.y;
  return is;
}

ostream &operator << (ostream &os,Point p){
  os<<fixed<<setprecision(12)<<p.x<<" "<<p.y;
  return os;
}


double area(Polygon s){
//求polygon面积 边给的都是连着的
  double res=0;
  for(int i=0;i<(int)s.size();i++){
    res+=Cross(s[i],s[(i+1)%s.size()])/2.0;
  }
  return res;
}

//判断点与多边形的位置
int contains(Polygon g,Point p){
  int n=g.size();
  bool x=false;
  for(int i=0;i<n;i++){
    Point a=g[i]-p,b=g[(i+1)%n]-p;
    if(fabs(Cross(a,b))<eps&& Dot(a,b)<eps) return 1;
    if(a.y>b.y) swap(a,b);
    if(a.y<eps&&eps<b.y&&Cross(a,b)>eps) x = !x;
  }
  return (x?2:0);
}

//判断凸包  所有点是逆时针顺序给出的
bool isConvex(Polygon p){
  bool f=1;
  int n=p.size();
  for(int i=0;i<n;i++){
    int t=ccw(p[(i+n-1)%n],p[i],p[(i+1)%n]);
    f&=t!=CLOCKWISE;
  }
  return f;
}

Polygon andrewScan(Polygon s){
    Polygon u,l;
    if(s.size()<3)return s;
    sort(s.begin(),s.end());
//这里注意是否需要去重
    u.push_back(s[0]);
    u.push_back(s[1]);
    l.push_back(s[s.size()-1]);
    l.push_back(s[s.size()-2]);
    for(int i=2;i<s.size();i++){
        for(int n=u.size();n>=2 && ccw(u[n-2],u[n-1],s[i])==1;n--){
            u.pop_back();
        }
        u.push_back(s[i]);
    }

    for(int i=s.size()-3;i>=0;i--){
        for(int n=l.size();n>=2 && ccw(l[n-2],l[n-1],s[i])==1;n--){
            l.pop_back();
        }
        l.push_back(s[i]);
    }
    reverse(l.begin(),l.end());
    for(int i=u.size()-2;i>=1;i--){
        l.push_back(u[i]);
    }
    return l;
}

double diameter(Polygon s){
//凸包直径   注意凸包都是逆时针的 然后直径是最长两点的长度
//也就是旋转卡壳
  Polygon p=s;
  int n=p.size();
  if(n==2) return Abs(p[0]-p[1]);
  int i=0,j=0;
  for(int k=0;k<n;k++){
    if(p[i]<p[k]) i=k;
    if(!(p[j]<p[k])) j=k;
  }
  //cout<<i<<" "<<j<<endl;
  double res=0;
  int si=i,sj=j;
  while(i!=sj||j!=si){
    //cout<<i<<" "<<j<<endl;
    res=max(res,Abs(p[i]-p[j]));
    if(Cross(p[(i+1)%n]-p[i],p[(j+1)%n]-p[j])<0.0){
      i=(i+1)%n;
    }else{
      j=(j+1)%n;
    }
  }
  return res;
}


Polygon convexCut(Polygon p,Line l){
//左半部分面积 可以调整至右半部分面积
  Polygon q;
  for(int i=0;i<(int)p.size();i++){
    Point a=p[i],b=p[(i+1)%p.size()];
    if(ccw(l.p1,l.p2,a)!=-1) q.push_back(a);
    if(ccw(l.p1,l.p2,a)*ccw(l.p1,l.p2,b)<0)
      q.push_back(getCrossPointLL(Line(a,b),l));
  }
  return q;
}

Point Polygon_center(Polygon p)
//求多边形重心
{
    int n=p.size();
    Point ans(0,0);
    if(area(p)==0) return ans;
    for(int i=0;i<n;i++)
    ans=ans+(p[i]+p[(i+1)%n])*Cross(p[i],p[(i+1)%n]);
    return ans/area(p)/6;
}
圆
struct Circle{
  Point c;
  double r;
  Circle(){}
  Circle(Point c,double r):c(c),r(r){}
};

两圆位置关系
int intersectCC(Circle c1,Circle c2){
  if(c1.r<c2.r) swap(c1,c2);
  double d=Abs(c1.c-c2.c);
  double r=c1.r+c2.r;
  if(equals(d,r)) return ICC_CIRCUMSCRIBE;
  if(d>r) return ICC_SEPERATE;
  if(equals(d+c2.r,c1.r)) return ICC_INSCRIBE;
  if(d+c2.r<c1.r) return ICC_CONTAIN;
  return ICC_INTERSECT;
}

pair<Point,Point> getCrossPoints(Circle c,Line l){
  Vector pr=project(c.c,l);
  Vector e=(l.p2-l.p1)/Abs(l.p2-l.p1);
  double base=sqrt(c.r*c.r-Norm(pr-c.c));
  return make_pair(pr+e*base,pr-e*base);
}//求圆和线的交点

double Arg(Vector p){return atan2(p.y,p.x);}
//atan2函数返回的是原点至点(x,y)的方位角，即与 x 轴的夹角。也可以理解为复数 x+yi 的辐角 [-pi,pi]
Vector polar(double a,double r){ return Point(cos(r)*a,sin(r)*a);}

pair<Point,Point> getCrossPoints(Circle c1,Circle c2){
  double d=Abs(c1.c-c2.c);
  double a=acos((c1.r*c1.r+d*d-c2.r*c2.r)/(2*c1.r*d));
  double t=Arg(c2.c-c1.c);
  return make_pair(c1.c+polar(c1.r,t+a),c1.c+polar(c1.r,t-a));
}//求圆和圆的交点（pair版本）

Vector rot(Vector a,double theta){//rotate.
    double x1 = cos(theta) * a.x - sin(theta) * a.y;
    double y1 = sin(theta) * a.x + cos(theta) * a.y;
    return Point(x1,y1);
}
vector <Point> tangent(Circle o,Point p){//求圆的切线，确保切线存在
    Vector a = p - o.c;
    double theta = acos(o.r / a.abs());
    Vector b = a / a.abs() * o.r;
    vector <Point> res;
    res.push_back(o.c + rot(b,theta));
    res.push_back(o.c + rot(b,-theta));
    return res;
}

//求在圆o1上的公切线
vector <Point> commonTangents(Circle o1,Circle o2){
    double rdif = o1.r - o2.r,rsum = o1.r + o2.r;
    Vector a = o2.c - o1.c,b = a / a.abs() * o1.r;
    vector <Point> res;
    res.clear();
    double alpha;
    if(abs(abs(rdif) - a.abs()) < eps){
        if(rdif < 0) alpha = acos(-1);
        else alpha = acos(1);
        res.push_back(o1.c + rot(b,alpha));
    }else if(abs(rdif) < a.abs() - eps){
        alpha = acos(rdif / a.abs());
        res.push_back(o1.c + rot(b,alpha));
        res.push_back(o1.c + rot(b,-alpha));
    }
    if(abs(rsum - a.abs()) < eps){
        alpha = acos(1);
        res.push_back(o1.c + rot(b,alpha));
    }else if(rsum < a.abs() - eps){
        alpha = acos(rsum / a.abs());
        res.push_back(o1.c + rot(b,alpha));
        res.push_back(o1.c + rot(b,-alpha));
    }
    return res;
}

//圆的中垂线 垂的是Line(p1,p2)
Line bisector(Point p1,Point p2){
  Circle c1=Circle(p1,Abs(p1-p2)),c2=Circle(p2,Abs(p1-p2));
  pair<Point,Point> p=getCrossPoints(c1,c2);
  if(Cross(p2-p1,p.first-p1)>0) swap(p.first,p.second);
  return Line(p.first,p.second);
}

/*
// 扇形面积，圆心在c，扇形a到b
double fan_area(double r, Point a, Point b){
    double angle = acos(Dot(a, b) / Abs(a) / Abs(b));
    if(sgn(Cross(a, b)) < 0) angle = -angle;
    return r * r * angle / 2;
}

// 求直线和圆交点
pair<Point, Point> intersection_of_line_and_circle (Circle c, Point l1, Point l2){
    Point p = project(c.c, Line(l1, l2));
    Point e = (l2 - l1) / Abs(l2 - l1);
    double base = sqrt(c.r * c.r - Abs(p - c.c) * Abs(p - c.c));
    return {p - e * base, p + e * base};
}

// 求三角形与圆相交面积，圆心(0,0)，三角形一定点过圆心
double intersection_area_of_triangle_and_circle (Point a, Point b, double r){
    Point p =project({0, 0},Line(a, b));
    double lena = Abs(a), lenb = Abs(b), lenp = Abs(p);
    if(!Point_on_seg(p,Line(a,b))) lenp = min(Abs(a), Abs(b));
    if(Dcmp(r, lena) >= 0 && Dcmp(r, lenb) >= 0) return Cross(a, b) / 2;
    if(sgn(Cross(a, b)) == 0) return 0;
    pair<Point, Point> P = intersection_of_line_and_circle(Circle({0, 0}, r), a, b);
    Point pa = P.first, pb = P.second;
    if(Dcmp(r, lenp) <= 0) return fan_area(r, a, b);
    if(Dcmp(r, lena) >= 0) return Cross(a, pb) / 2 + fan_area(r, pb, b);
    if(Dcmp(r, lenb) >= 0) return Cross(pa, b) / 2 + fan_area(r, a, pa);
    return fan_area(r, a, pa) + Cross(pa, pb) / 2 + fan_area(r, pb, b);
}

// 多边形与圆相交面积
double intersection_area_of_polygon_and_circle(Polygon p, Circle c){
    double res = 0;
    int n=p.size();
    for(int i = 0; i < n; i++){
        double rest = intersection_area_of_triangle_and_circle(p[i], p[(i+1) % n], c.r);
        res += rest;
    }
    return res;
}
*/凸多边形和圆交大板块


//两圆相交面积
double intersection_area_of_two_circles(Circle c1, Circle c2){
    double d = Abs(c1.c - c2.c);
    if(c1.r + c2.r < d + eps) return 0;
    if(d < fabs(c1.r - c2.r) + eps) {
        double r = min(c1.r, c2.r);
        return pi * r * r;
    }
    double a1 = acos((d * d + c1.r * c1.r - c2.r * c2.r) / (2 * c1.r * d));
    double a2 = acos((d * d + c2.r * c2.r - c1.r * c1.r) / (2 * c2.r * d));
    double s1 = c1.r * c1.r * sin(2 * a1) / 2 + c2.r * c2.r * sin(2 * a2) / 2;
    double s2 = c1.r * c1.r * a1 + c2.r * c2.r * a2;
    return s2 - s1;
}

//三角形内接圆
Circle inscribed_circle_of_triangle(Point p1, Point p2, Point p3){
    Point o; double r;
    double a = Abs(p1-p2), b = Abs(p2-p3), c = Abs(p1-p3);
    double area = fabs(Cross(p3 - p1, p3 - p2))/ 2.0;
    r = 2 * area / (a + b + c);
    o.x = (a * p3.x + b * p1.x + c * p2.x) / (a + b + c);
    o.y = (a * p3.y + b * p1.y + c * p2.y) / (a + b + c);
    return Circle(o, r);
}

// 求三角形外接圆
Circle circumcircle_of_triangle(Point p1, Point p2, Point p3) {  //三点定圆
    Point o; double r;
    double a, b, c, d, e, f;
    a = p1.x - p2.x, b = p1.y - p2.y, c = p1.x - p3.x, d = p1.y - p3.y;
    e = (p1.x * p1.x - p2.x * p2.x) - (p2.y * p2.y - p1.y * p1.y);
    f = (p1.x * p1.x - p3.x * p3.x) - (p3.y * p3.y - p1.y * p1.y);
    o.x = (b * f - d * e) / (2 * b * c - 2 * a * d);
    o.y = (c * e - a * f) / (2 * b * c - 2 * a * d);
    r = Abs(o-p1);
    return Circle(o, r);
}

极角排序
//从第三象限开始
    vector<Vec> a(n);
    for(int i=0;i<n;i++)
        scanf("%lld%lld", &a[i].x, &a[i].y);
    auto quad = [](Vec a) {
        if (a.y < 0) return 1;
        if (a.y > 0) return 4;
        if (a.x < 0) return 5;
        if (a.x > 0) return 3;
        return 2;
    };
    sort(a.begin(), a.end(), [&](Vec a, Vec b) {
        int qa = quad(a), qb = quad(b);
        return qa == qb ? a.cross(b) > 0 : qa < qb;
    });
    //  for(int i=0;i<n;i++)
    // printf("%lld %lld \n", a[i].x, a[i].y);

//从第一象限开始
LL cross(Point a, Point b){
    return (LL)a.x * b.y - (LL)a.y * b.x;
}

bool up(Point a){
    return a.y > 0 || (a.y == 0 && a.x >= 0);
}
  sort(p + 1, p + n + 1, [&](Point a, Point b){
            if (up(a) != up(b)) return up(a) > up(b);
            return cross(a, b) > 0;
        });


半平面交
Vector Normal(Vector A){   double L = Abs(A);  return Vector(-A.y / L, A.x / L);}

//有向直线。它的左边就是对应的半平面
struct Line
{
    Point P;//直线上任意一点
    Vector v;//方向向量。左边就是对应的半平面
    double deg;//极角
    Line(){}
    Line(Point P, Vector v):P(P), v(v){deg = atan2(v.y, v.x);}
    bool operator < (const Line& L)const {//排序时使用的比较运算符
        return deg < L.deg;
    }
};
//半平面交
//半平面交一般是一个凸多边形，但是有时候会是一个无界多边形
//甚至会是一个直线、线段、点，但是结果一定是凸的

//点p在有向直线L的左边（线上的不算）（叉积大于0a在b的左侧，小于0在右侧[sin夹角]）
bool on_left(Line L, Point P){return Cross(L.v, P - L.P) > 0;}
//两个有向直线的交点/假定交点唯一存在
Point get_intersection(Line a, Line b)
{
    Vector u = a.P - b.P;
    double t = Cross(b.v, u) / Cross(a.v, b.v);
    return a.P + a.v * t;
}
//半平面交的主过程
int half_plane_intersection(Line* L, int n, Point* poly)
{
    sort(L, L + n);//按照极角排序

    int first, last;//双端队列
    Point *p = new Point[n];//p[i]为q[i]和q[i + 1]的交点
    Line *q = new Line[n];//手写的Line类型的双端队列（数组）
    q[first = last = 0] = L[0];//双端队列初始化的时候只有一个半平面L[0]
    for(int i = 1; i < n; ++ i){
        while(first < last && !on_left(L[i], p[last - 1]))last -- ;
        while(first < last && !on_left(L[i], p[first]))first ++ ;
        q[++ last] = L[i];//新的点是一定要放进去的
        if(fabs(Cross(q[last].v, q[last - 1].v)) < eps){
        //相邻的两个向量平行且同向，取内侧的那一个
            last -- ;//如果新的向量上的一个点在老的向量的左侧就取新的
            if(on_left(q[last], L[i].P))q[last] = L[i];
        }
        if(first < last)p[last - 1] = get_intersection(q[last - 1], q[last]);
    }
    while(first < last && !on_left(q[first], p[last - 1]))last -- ;
    //删除无用的平面

    if(last - first <= 1)return 0;//空集
    p[last] = get_intersection(q[last], q[first]);//计算首尾两个半平面交（环状）
    //从手写deque中复制答案到输出数组中
    int m = 0;
    for(int i = first ; i <= last; ++ i)poly[m ++ ] = p[i];
    return m;
}



闵科夫斯基和
Vector V1[N],V2[N];
int Mincowski(Polygon P1,Polygon P2,Vector *V){//【闵可夫斯基和】求两个凸包{P1},{P2}的向量集合{V}={P1+P2}构成的凸包
    int n=P1.size();
    int m=P2.size();
    for(int i=1;i<=n;++i)V1[i]=P1[i<n?i+1:1]-P1[i];
    for(int i=1;i<=m;++i)V2[i]=P2[i<m?i+1:1]-P2[i];
    int t=0,i=1,j=1;V[++t]=P1[1]+P2[1];
    while(i<=n&&j<=m)++t,V[t]=V[t-1]+(Dcmp(Cross(V1[i],V2[j]))>0?V1[i++]:V2[j++]);
    while(i<=n)++t,V[t]=V[t-1]+V1[i++];
    while(j<=m)++t,V[t]=V[t-1]+V2[j++];
    return t;
}


多边形内网格点个数
#define abs(x) ((x)>0?(x):-(x))
struct point{int x,y;};
int gcd(int a,int b){return b?gcd(b,a%b):a;}
//多边形上的网格点个数
int grid_onedge(int n,point* p){
    int i,ret=0;
    for (i=0;i<n;i++)
        ret+=gcd(abs(p[i].x-p[(i+1)%n].x),abs(p[i].y-p[(i+1)%n].y));
    return ret;
}
//多边形内的网格点个数
int grid_inside(int n,point* p){
    int i,ret=0;
    for (i=0;i<n;i++)
        ret+=p[(i+1)%n].y*(p[i].x-p[(i+2)%n].x);
    return (abs(ret)-grid_onedge(n,p))/2+1;
}


多边形并面积 板子
#include<bits/stdc++.h>
using namespace std;
#define mp make_pair
typedef long long ll;
const double inf=1e200;
const double eps=1e-12;
const double pi=4*atan(1.0);
int dcmp(double x){ return fabs(x)<eps?0:(x<0?-1:1);}
struct point{
    double x,y;
    point(double a=0,double b=0):x(a),y(b){}
};
point operator +(point A,point B) { return point(A.x+B.x,A.y+B.y);}
point operator -(point A,point B) { return point(A.x-B.x,A.y-B.y);}
point operator *(point A,double p){ return point(A.x*p,A.y*p);}
point operator /(point A,double p){ return point(A.x/p,A.y/p);}
bool operator ==(const point& a,const point& b){
    return fabs(a.x-b.x)<eps&&fabs(a.y-b.y)<eps;
}
double dot(point A,point B){ return A.x*B.x+A.y*B.y;}
double det(point A,point B){ return A.x*B.y-A.y*B.x;}
double det(point O,point A,point B){ return det(A-O,B-O);}
double length(point A){ return sqrt(dot(A,A));}
double area(vector<point>p){
    double ans=0; int sz=p.size();
    for(int i=1;i<sz-1;i++) ans+=det(p[i]-p[0],p[i+1]-p[0]);
    return ans/2.0;
}
double seg(point O,point A,point B){
    if(dcmp(B.x-A.x)==0) return (O.y-A.y)/(B.y-A.y);
    return (O.x-A.x)/(B.x-A.x);
}
vector<point>pp[110];
pair<double,int>s[110*60];
double polyunion(vector<point>*p,int N){
    double res=0;
    for(int i=0;i<N;i++){
        int sz=p[i].size();
        for(int j=0;j<sz;j++){
            int m=0;
            s[m++]=mp(0,0);
            s[m++]=mp(1,0);
            point a=p[i][j],b=p[i][(j+1)%sz];
            for(int k=0;k<N;k++){
                if(i!=k){
                    int sz2=p[k].size();
                    for(int ii=0;ii<sz2;ii++){
                        point c=p[k][ii],d=p[k][(ii+1)%sz2];
                        int c1=dcmp(det(b-a,c-a));
                        int c2=dcmp(det(b-a,d-a));
                        if(c1==0&&c2==0){
                            if(dcmp(dot(b-a,d-c))){
                                s[m++]=mp(seg(c,a,b),1);
                                s[m++]=mp(seg(c,a,b),-1);
                            }
                        }
                        else{
                            double s1=det(d-c,a-c);
                            double s2=det(d-c,b-c);
                            if(c1>=0&&c2<0) s[m++]=mp(s1/(s1-s2),1);
                            else if(c1<0&&c2>=0) s[m++]=mp(s1/(s1-s2),-1);
                        }
                    }
                }    
            }
            sort(s,s+m);
            double pre=min(max(s[0].first,0.0),1.0),now,sum=0;
            int cov=s[0].second;
            for(int j=1;j<m;j++){
                now=min(max(s[j].first,0.0),1.0);
                if(!cov) sum+=now-pre;
                cov+=s[j].second;
                pre=now;
            }
            res+=det(a,b)*sum;
        }
    }
    return res/2;
}
int main()
{
    int N,M,i,j; point tp;
    scanf("%d",&N);
    for(i=0;i<N;i++){
        scanf("%d",&M);
        for(j=0;j<M;j++){
            scanf("%lf%lf",&tp.x,&tp.y);
            pp[i].push_back(tp);
        }
    }
    double t1=0,t2=polyunion(pp,N);
    for(i=0;i<N;i++) t1+=area(pp[i]);
    printf("%.7lf %.7lf\n",-t1,-t2);
    return 0;
}
