#include <bits/stdc++.h>
using namespace std;
double INF =  1e100;
double EPS = 1e-12;

struct PT{
    double x, y;
    PT() {}
    PT(double x, double y) : x(x), y(y) {}
    PT(const PT &p) : x(p.x), y(p.y) {}
    PT operator + (const PT &p) const {return PT(x+p.x, y+p.y); }
    PT operator - (const PT &p) const {return PT(x-p.x, y-p.y); }
    PT operator * (double c) const {return PT(x*c, y*c); }
    PT operator / (double c) const {return PT(x/c, y/c); }
    double dist(){
        return sqrt(x*x + y*y);
    }
    double dist(const PT &p){
        double xx = x - p.x, yy = y - p.y;
        return sqrt(xx*xx + yy*yy);
    }
};

double dot(PT p, PT q) {return p.x*q.x + p.y*q.y;}
double dist2(PT p, PT q) {return dot(p-q,p-q); }
double cross(PT p, PT q) {return p.x*q.y - p.y*q.x;}

// CCW check
double CCW(PT a, PT b, PT c){ return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);}
// rotate a point
PT RotateCCW90(PT p) {return PT(-p.y, p.x); }
PT RotateCW90(PT p) {return PT(p.y, -p.x); }
PT RotateCCW(PT p, double t){
    return PT(p.x*cos(t)-p.y*sin(t) , p.x*sin(t)+p.y*cos(t) );
}

// project point c onto line through a and b
// assuming a!= b
PT ProjectPointLine(PT a, PT b, PT c){
    return a + (b-a)*dot(c-a,b-a)/dot(b-a,b-a);
}

// project point c onto line segment through a and b
PT ProjectPointSegment(PT a, PT b, PT c){
    double r = dot(b-a,b-a);
    if(fabs(r) < EPS) return a;
    r = dot(c-a,b-a)/r;
    if(r < 0) return a;
    if(r > 1) return b;
    return a + (b-a)*r;
}

// compute ditance from c to segment between a and b
double DistancePointSegment(PT a, PT b, PT c){
    return sqrt(dist2(c,ProjectPointSegment(a,b,c)));
}

// compute distance between point (x,y,x) and plane ax+by+cy= d
double DistancePointPlane(double x, double y, double z, double a, double b, double c, double d){
    return fabs(a*x+b*y+c*z-d)/sqrt(a*a+b*b+c*c);
}

// determine if lines from a to b and c to d are parallel or collinear
bool LinesParallel(PT a, PT b, PT c, PT d) {
    return fabs(cross(b-a, c-d)) < EPS;
}

bool LinesCollinear(PT a, PT b, PT c, PT d) {
    return LinesParallel(a, b, c, d) && fabs(cross(a-b, a-c)) < EPS && fabs(cross(c-d, c-a)) < EPS;
}

// determine if line segment from a to b and from c to d
bool SegmentsIntersect(PT a, PT b, PT c, PT d){
    if(LinesCollinear(a,b,c,d)){
        if(dist2(a,c) < EPS || dist2(a,b) < EPS || dist2(b,c) < EPS || dist2(b,d) < EPS) return true;
        if(dot(c-a,c-b) > 0 && dot(d-a,d-b) > 0 && dot(c-b,d-b) > 0) return false;
        return true;
    }
    if(cross(d-a,b-a) * cross(c-a, b-a) > 0) return false;
    if(cross(a-c,d-c) * cross(b-c, d-c) > 0) return false;
    return true;
}

// compute intersection of line passing through a and b
// with line passing through c and d, assuming that unique
// intersection exists; for segment intersection, check if
// segments intersect first
PT ComputeLineIntersection(PT a, PT b, PT c, PT d) {
    b=b-a; d=c-d; c=c-a;
    assert(dot(b, b) > EPS && dot(d, d) > EPS);
    return a + b*cross(c, d)/cross(b, d);
}

// compute outter center of circle given three points
PT OutterCircleCenter1(PT a, PT b, PT c){
    b = (a+b)/2;
    c = (a+c)/2;
    return ComputeLineIntersection(b,b+RotateCW90(a-b),c, c+RotateCW90(a-c));
}

PT OutterCircleCenter2(PT a, PT b, PT c){
    double dx1, dx2, dy1, dy2, dx3, dy3;
    dx1 = b.x - a.x;
    dx2 = c.x - a.x;
    dy1 = b.y - a.y;
    dy2 = c.y - a.y;
    dx3 = c.x - b.x;
    dy3 = c.y - b.y;
    double tq = 2 * (dy1*dx2 - dx1*dy2);
    double tp = -dx2*dx3 - dy2*dy3;
    if (abs(tq) < EPS) {
        return PT(-1, -1);
    }
    double sq = tq;
    double sp = (a.x + b.x) * tq / 2 - (dy1)*tp;
    double sr = (a.y + b.y) * tq / 2 + (dx1)*tp;
    return PT(sp / sq, sr / sq);
}

// determine if point is in a possibly non-convex polygon (by William
// Randolph Franklin); returns 1 for strictly interior points, 0 for
// strictly exterior points, and 0 or 1 for the remaining points.
// Note that it is possible to convert this into an *exact* test using
// integer arithmetic by taking care of the division appropriately
// (making sure to deal with signs properly) and then by writing exact
// tests for checking point on polygon boundary
bool PointInPolygon(const vector<PT>&p, PT q) {
    bool c = 0;
    for (int i = 0; i <p.size(); i++) {
        int j = (i + 1) % p.size();
        if (( (p[i].y <= q.y &&q.y <p[j].y) || (p[j].y <= q.y &&q.y <p[i].y)) && q.x <p[i].x + (p[j].x - p[i].x) * (q.y - p[i].y) / (p[j].y - p[i].y))
            c = !c;
    }
    return c;
}

// determine if point is on the boundary of a polygon
bool PointOnPolygon(const vector<PT>&p, PT q) {
    for (int i = 0; i <p.size(); i++)
        if (dist2(ProjectPointSegment(p[i], p[(i + 1) % p.size()], q), q) < EPS)
            return true;
    return false;
}

// compute intersection of line through points a and b with
// circle centered at c with radius r > 0
vector<PT> CircleLineIntersection(PT a, PT b, PT c, double r) {
    vector<PT> ret;
    b=b-a;
    a=a-c;
    double A = dot(b, b);
    double B = dot(a, b);
    double C = dot(a, a) - r*r;
    double D = B*B - A*C;
    if (D < -EPS) return ret;
    ret.push_back(c+a+b*(-B + sqrt(D + EPS)) / A);
    if (D > EPS)
        ret.push_back(c+a+b*(-B - sqrt(D)) / A);
    return ret;
}

// compute intersection of circle centered at a with radius r
// with circle centered at b with radius R
vector<PT> CircleCircleIntersection(PT a, PT b, double r, double R) {
    vector<PT> ret;
    double d = sqrt(dist2(a, b));
    if (d >r + R || d + min(r, R) < max(r, R)) return ret;
    double x = (d*d - R*R + r*r) / (2 * d);
    double y = sqrt(r*r - x*x);
    PT v = (b-a) / d;
    ret.push_back(a+ v*x + RotateCCW90(v)*y);
    if (y > 0)
        ret.push_back(a+ v*x - RotateCCW90(v)*y);
    return ret;
}

//circle b,R and tangencial from a
vector<PT> CircleTangencial(PT a, PT b, double r)
{
    PT res;
    double dd = a.dist(b);
    double ang = acos(r / dd);
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double zz;
    vector<PT> rr;
    res.x = dx * cos(ang) - dy * sin(ang);
    res.y = dx * sin(ang) + dy * cos(ang);
    zz = res.dist();
    res.x = res.x / zz * r;
    res.y = res.y / zz * r;
    res= res +b;
    rr.push_back(res);
    ang = -ang;
    res.x = dx * cos(ang) - dy * sin(ang);
    res.y = dx * sin(ang) + dy * cos(ang);
    zz = res.dist();
    res.x = res.x / zz * r;
    res.y = res.y / zz * r;
    res= res +b;
    rr.push_back(res);
    return rr;
}


// This code computes the area or centroid of a (possibly nonconvex)
// polygon, assuming that the coordinates are listed in a clockwise or
// counterclockwise fashion. Note that the centroid is often known as
// the "center of gravity" or "center of mass".
double ComputeSignedArea(const vector<PT> &p) {
    double area = 0;
    for(int i = 0; i < p.size(); i++) {
        int j = (i+1) % p.size();
        area += p[i].x*p[j].y - p[j].x*p[i].y;
    }
    return area / 2.0;
}

double ComputeArea(const vector<PT> &p) {
    return fabs(ComputeSignedArea(p));
}

PT ComputeCentroid(const vector<PT> &p) {
    PT c(0,0);
    double scale = 6.0 * ComputeSignedArea(p);
    for (int i = 0; i < p.size(); i++){
        int j = (i+1) % p.size();
        c = c + (p[i]+p[j])*(p[i].x*p[j].y - p[j].x*p[i].y);
    }
    return c / scale;
}

// tests whether or not a given polygon (in CW or CCW order) is simple
bool IsSimple(const vector<PT> &p) {
    for (int i = 0; i < p.size(); i++) {
        for (int k = i+1; k < p.size(); k++) {
            int j = (i+1) % p.size();
            int l = (k+1) % p.size();
            if (i == l || j == k) continue;
            if (SegmentsIntersect(p[i], p[j], p[k], p[l])) return false;
        }
    }
    return true;
}

int main(){

}
