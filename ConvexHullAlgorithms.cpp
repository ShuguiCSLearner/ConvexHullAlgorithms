// ConvexHullAlgorithms.cpp : Defines the entry point for the application.
//

#include "framework.h"
#include "ConvexHullAlgorithms.h"
#include <string>

#include <cmath>
#include <stack> // for Graham Scan Code Hull Algorithm
#include <list>
#include <iostream>
#include <stdlib.h> // qsort
#include <cstdlib> // random number generator
#include <ctime> // for random number
using namespace std;

#define POINT_CONVEX_HULL 1
#define QUICKHULL 2
#define MINKOWSKI_DIFFERENCE 3
#define MINKOWSKI_SUM 4
#define GJK 5


POINT convexHullPoints[13];
// Size of convexHullList
int listSize = 0;

// array that holds random vertices
int numberOfVertex = 12;
POINT randomVertices[12];

// holds the reference of one vertex in the convex hull list
POINT v0;

// center vertex reference
POINT centerC;


// ------------------------------------------------minkowski sum global data----------------------------------------------//
// when computing minkoski sum 
// call -> MinkowskiSumGrahamScanCodeForHullA(minkowskiSumUnsortedVerticesA, 6);   for computing convex hull A
// then -> MinkowskiSumGrahamScanCodeForHullB(minkowskiSumUnsortedVerticesB, 6);   for computing convex hull B
// then -> MinkowskiSum(x1,x2,y1,y1) where x1 is the left border, x2 is the right border, y1 is the lower border, and y2 is the upport border of your window.
// then -> MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);   compute the third convex hull C = A + B.
// then -> poly()
// For GrahanScanCode Function
// above 3 point dont care since used by graham scan code, which can be replace by local variable as well. Not a issue but I want to use global...
POINT minkowskiSumUnsortedVerticesA[6];
POINT minkowskiSumUnsortedVerticesB[6];
POINT minkowskiSumUnsortedVerticesC[36];
// above is the unsorted vertices
// below is sorted convex hull points
POINT minkowskiSumConvexHullPointsA[7];
POINT minkowskiSumConvexHullPointsB[7];
POINT minkowskiSumConvexHullPointsC[37];
int minkowskiSumConvexHullListSizeA = 0;
int minkowskiSumConvexHullListSizeB = 0;
int minkowskiSumConvexHullListSizeC = 0;
// ----------------------------------------------minkowski sum global data end--------------------------------------------//

// return the second vertex of the stack
POINT SecondVertexInStack(stack<POINT>& S) {
    POINT first = S.top();
    S.pop();
    POINT second = S.top();
    S.push(first); // push the popped vertex back
    return second;
}

// swap 2 point by reference
void Swap(POINT& v1, POINT& v2) {
    POINT temp = v1;
    v1 = v2;
    v2 = temp;
}

// return the square of distance between v1 and v2: positive number // hypotenuse
double DistanceSquare(POINT v1, POINT v2) {
    return ((double)v1.x - (double)v2.x) * ((double)v1.x - (double)v2.x) + ((double)v1.y - (double)v2.y) * ((double)v1.y - (double)v2.y);
}

// check whether vertex c is on the clockwise or counterclockwise or collinear with segment ab
// 0 = collinear, -1 = ccw, 1 = cw
int Direction(POINT a, POINT b, POINT c) {
    double area = ((double)b.y - (double)a.y) * ((double)c.x - (double)b.x) - ((double)b.x - (double)a.x) * ((double)c.y - (double)b.y);
    if (area == 0) return 0;
    if (area > 0) return 1; // cw
    else return -1; // ccw
}

// compare whether v2 is in clockwise of segment v0v1 or ccw or collinear
int Compare(const void* vp1, const void* vp2) {
    // type casting
    POINT* v1 = (POINT*)vp1;
    POINT* v2 = (POINT*)vp2;

    int d = Direction(v0, *v1, *v2);
    if (d == 0) { // case of collinear
        if (DistanceSquare(v0, *v2) >= DistanceSquare(v0, *v1)) {
            return -1;
        }
        else return 1;
    }
    return d; // -1 is ccw and 1 = cw
}

// compute the convex hull
void GrahamScanCode(POINT vertices[], int n) {
    for (int i = 0; i < 12; i++) {
        convexHullPoints[i].x = vertices[i].x;
        convexHullPoints[i].y = vertices[i].y;
    }
    double minY = convexHullPoints[0].y;
    int minIndex = 0;
    // find bottommost left vertex
    for (int i = 1; i < n; i++) {
        double currentY = convexHullPoints[i].y;
        if ((currentY < minY) || (minY == currentY && convexHullPoints[i].x < convexHullPoints[minIndex].x)) {
            minY = currentY;
            minIndex = i;
        }
    }
    // the first vertex is the bottommost left vertex
    Swap(convexHullPoints[0], convexHullPoints[minIndex]);
    v0 = convexHullPoints[0];
    qsort(&convexHullPoints[1], n - 1, sizeof(POINT), Compare);

    int j = 1;
    for (int i = 1; i < n; i++) {
        while (i < n - 1 && Direction(v0, convexHullPoints[i], convexHullPoints[i + 1]) == 0) {
            i++;
        }
        convexHullPoints[j] = convexHullPoints[i];
        j++;
    }
    if (j < 3) return;
    stack<POINT> vertexStack;
    vertexStack.push(convexHullPoints[0]);
    vertexStack.push(convexHullPoints[1]);
    vertexStack.push(convexHullPoints[2]);
    for (int i = 3; i < j; i++) {
        while (vertexStack.size() > 1 && Direction(SecondVertexInStack(vertexStack), vertexStack.top(), convexHullPoints[i]) != -1) {
            vertexStack.pop();
        }
        vertexStack.push(convexHullPoints[i]);
    }
    // empty the convexHullList in global if it contain elements
    while (listSize != 0) {
        listSize--;
    }
    int counter = 0;
    while (!vertexStack.empty()) {
        POINT v = vertexStack.top();
        convexHullPoints[counter] = v;
        counter++;
        //convexHullList.push_back(v);
        listSize++;
        //cout << "(" << v.x << ", " << v.y << ")" << endl;
        vertexStack.pop();
    }
    convexHullPoints[counter] = convexHullPoints[0];
}


// from Joseph O'Rourke's computational geometry
// finding  and return the area for the triangle form by 3 vertices.
double Area2(POINT a, POINT b, POINT c) {
    return ((double)b.x - (double)a.x) * ((double)c.y - (double)a.y) - ((double)c.x - (double)a.x) * ((double)b.y - (double)a.y);
}

// check if vertex a,b,c are collinear
bool IsCollinear(POINT a, POINT b, POINT c) {
    return (Area2(a, b, c) == 0);
}

// check if vertex c is between vertex a and b.
bool IsBetween(POINT a, POINT b, POINT c) {
    if (!IsCollinear(a, b, c)) return false;
    if ((double)a.x != (double)b.x) { // vertex a and b are not horizontal line
        return ((((double)a.x <= (double)c.x) && ((double)c.x <= (double)b.x)) || (((double)a.x >= (double)c.x) && ((double)c.x >= (double)b.x)));
    }
    else { // vertex a and b form a horizontal line
        return ((((double)a.y <= (double)c.y) && ((double)c.y <= (double)b.y)) || (((double)a.y >= (double)c.y) && ((double)c.y >= (double)b.y)));
    }
}

// find the larger number
double Max(double a, double b) {
    if (a > b) return a;
    else return b;
}

// find the smaller number
double Min(double a, double b) {
    if (a < b) return a;
    else return b;
}

// check whether vertex c stretch horizontally toward right that form a single point segment intersect with segment ab or not.
bool IsIntersect(POINT a, POINT b, POINT c) {
    if (IsBetween(a, b, c)) return false; // vertex c lays on the edge of the polygon
    if ((double)c.x < Max((double)a.x, (double)b.x)) {
        if ((double)c.y >= Min((double)a.y, (double)b.y) && c.y < Max((double)a.y, (double)b.y)) {  // y coordinate of c within the range of segment ab's y coordinates
            if ((double)a.y != (double)b.y) {
                // (x2 - x1) / (y2 - y1) * (y3 - y1) + x1
                double xintersection = ((double)b.x - (double)a.x) / ((double)b.y - (double)a.y) * ((double)c.y - (double)a.y) + (double)a.x;
                if ((double)c.x < xintersection) {
                    return true;
                }
            }
        }
    }
    return false;
}

// Vertex hull[] consists a list of vertex of the convex hull in order.
// Vertex c is the point we want to check if it it in the 
bool IsVertexInsideConvexHull(POINT hull[], POINT c, int hullSize) {
    int counter = 0;
    for (int i = 0; i < hullSize; i++) {
        //advance(it1, i);
        //advance(it2, (i + 1) % hullSize);
        POINT a = hull[i];
        POINT b = hull[(i + 1) % hullSize];
        if (IsBetween(a, b, c)) return false; // c is on edge of segment ab
        cout << "i: " << i << endl;
        cout << "a: " << a.x << "," << a.y << endl;
        cout << "b: " << b.x << "," << b.y << endl;
        if (IsIntersect(a, b, c)) {
            counter++;
        }
        cout << "current counter: " << counter << endl;
    }
    cout << "Counter: " << counter << endl;
    if (counter % 2 == 1) return true;
    else return false;
}

// This function will return a random double within range [min, max] inclusively
double GetRandomDouble(double min, double max) {
    int number = rand() % (int)(max - min + 1) + (int)min;
    return (double)number;
}

// Check if randomvertices[] contains vertex with same x and y up to n th index
bool ContainsVertex(double x, double y, int n) {
    for (int i = 0; i < n; i++) {
        if (randomVertices[i].x == x && randomVertices[i].y == y) return true;
    }
    return false;
}

// generating 12 random vertex
void GeneratingVertex(double xmin, double xmax, double ymin, double ymax) {
    for (int i = 0; i < numberOfVertex; i++) { // x = [400, 1300], y = [100, 700]
        double x = GetRandomDouble(xmin, xmax);
        double y = GetRandomDouble(ymin, ymax);
        while (ContainsVertex(x, y, i)) { // generate and no repeat vertex
            x = GetRandomDouble(xmin, xmax);
            y = GetRandomDouble(ymin, ymax);
        }
        randomVertices[i].x = x;
        randomVertices[i].y = y;
    }
}

/*

    ----------------------------------------------------- Minkowski Sum Algorithm's function -----------------------------------------------------------------

*/
// ******using minkowskiSumUnsortedVerticesA   -> vertices;   **int n = 6
// compute the convex hull for minkowskiSumHull A       // **int n = 6     minkowskiSumConvexHullListSizeA change = number of convex hull points
void MinkowskiSumGrahamScanCodeForHullA(POINT vertices[], int n) {
    for (int i = 0; i < n; i++) {
        minkowskiSumConvexHullPointsA[i].x = vertices[i].x;
        minkowskiSumConvexHullPointsA[i].y = vertices[i].y;
    }
    double minY = minkowskiSumConvexHullPointsA[0].y;
    int minIndex = 0;
    // find bottommost left vertex
    for (int i = 1; i < n; i++) {
        double currentY = minkowskiSumConvexHullPointsA[i].y;
        if ((currentY < minY) || (minY == currentY && minkowskiSumConvexHullPointsA[i].x < minkowskiSumConvexHullPointsA[minIndex].x)) {
            minY = currentY;
            minIndex = i;
        }
    }
    // the first vertex is the bottommost left vertex
    Swap(minkowskiSumConvexHullPointsA[0], minkowskiSumConvexHullPointsA[minIndex]);
    v0 = minkowskiSumConvexHullPointsA[0];                                        //     <- need a new global variable  POINT type   was named vo
    qsort(&minkowskiSumConvexHullPointsA[1], n - 1, sizeof(POINT), Compare);

    int j = 1;
    for (int i = 1; i < n; i++) {
        while (i < n - 1 && Direction(v0, minkowskiSumConvexHullPointsA[i], minkowskiSumConvexHullPointsA[i + 1]) == 0) {   // <- change v0
            i++;
        }
        minkowskiSumConvexHullPointsA[j] = minkowskiSumConvexHullPointsA[i];
        j++;
    }
    if (j < 3) return;
    stack<POINT> vertexStack;
    vertexStack.push(minkowskiSumConvexHullPointsA[0]);
    vertexStack.push(minkowskiSumConvexHullPointsA[1]);
    vertexStack.push(minkowskiSumConvexHullPointsA[2]);
    for (int i = 3; i < j; i++) {
        while (vertexStack.size() > 1 && Direction(SecondVertexInStack(vertexStack), vertexStack.top(), minkowskiSumConvexHullPointsA[i]) != -1) {
            vertexStack.pop();
        }
        vertexStack.push(minkowskiSumConvexHullPointsA[i]);
    }
    // empty the convexHullList in global if it contain elements
    while (minkowskiSumConvexHullListSizeA != 0) {                              //      <- need a new global variable  int type was named  listSize
        minkowskiSumConvexHullListSizeA--;                                      //      < change listSize
    }
    int counter = 0;
    while (!vertexStack.empty()) {
        POINT v = vertexStack.top();
        minkowskiSumConvexHullPointsA[counter] = v;                   //      <- need a new global variable  POINT[] type was name convexHullPoints
        counter++;
        minkowskiSumConvexHullListSizeA++;                                      //      <- change listSize
        vertexStack.pop();
    }
    minkowskiSumConvexHullPointsA[counter] = minkowskiSumConvexHullPointsA[0];     //      <- change convexHullPoints
}

// ******using minkowskiSumUnsortedVerticesB   -> vertices;   **int n = 6
// compute the convex hull for minkowskiSumHull B       // **int n = 6     minkowskiSumConvexHullListSizeB change = number of convex hull points
void MinkowskiSumGrahamScanCodeForHullB(POINT vertices[], int n) {
    for (int i = 0; i < n; i++) {
        minkowskiSumConvexHullPointsB[i].x = vertices[i].x;
        minkowskiSumConvexHullPointsB[i].y = vertices[i].y;
    }
    double minY = minkowskiSumConvexHullPointsB[0].y;
    int minIndex = 0;
    // find bottommost left vertex
    for (int i = 1; i < n; i++) {
        double currentY = minkowskiSumConvexHullPointsB[i].y;
        if ((currentY < minY) || (minY == currentY && minkowskiSumConvexHullPointsB[i].x < minkowskiSumConvexHullPointsB[minIndex].x)) {
            minY = currentY;
            minIndex = i;
        }
    }
    // the first vertex is the bottommost left vertex
    Swap(minkowskiSumConvexHullPointsB[0], minkowskiSumConvexHullPointsB[minIndex]);
    v0 = minkowskiSumConvexHullPointsB[0];                                        //     <- need a new global variable  POINT type   was named minkowskiSumHullBv0
    qsort(&minkowskiSumConvexHullPointsB[1], n - 1, sizeof(POINT), Compare);

    int j = 1;
    for (int i = 1; i < n; i++) {
        while (i < n - 1 && Direction(v0, minkowskiSumConvexHullPointsB[i], minkowskiSumConvexHullPointsB[i + 1]) == 0) {   // <- change minkowskiSumHullBv0
            i++;
        }
        minkowskiSumConvexHullPointsB[j] = minkowskiSumConvexHullPointsB[i];
        j++;
    }
    if (j < 3) return;
    stack<POINT> vertexStack;
    vertexStack.push(minkowskiSumConvexHullPointsB[0]);
    vertexStack.push(minkowskiSumConvexHullPointsB[1]);
    vertexStack.push(minkowskiSumConvexHullPointsB[2]);
    for (int i = 3; i < j; i++) {
        while (vertexStack.size() > 1 && Direction(SecondVertexInStack(vertexStack), vertexStack.top(), minkowskiSumConvexHullPointsB[i]) != -1) {
            vertexStack.pop();
        }
        vertexStack.push(minkowskiSumConvexHullPointsB[i]);
    }
    // empty the convexHullList in global if it contain elements
    while (minkowskiSumConvexHullListSizeB != 0) {                              //      <- need a new global variable  int type was named  minkowskiSumConvexHullListSizeB
        minkowskiSumConvexHullListSizeB--;                                      //      < change minkowskiSumConvexHullListSizeB
    }
    int counter = 0;
    while (!vertexStack.empty()) {
        POINT v = vertexStack.top();
        minkowskiSumConvexHullPointsB[counter] = v;                   //      <- need a new global variable  POINT[] type was name minkowskiSumConvexHullPointsB
        counter++;
        minkowskiSumConvexHullListSizeB++;                                      //      <- change minkowskiSumConvexHullListSizeB
        vertexStack.pop();
    }
    minkowskiSumConvexHullPointsB[counter] = minkowskiSumConvexHullPointsB[0];     //      <- change minkowskiSumConvexHullPointsB
}

// ******using minkowskiSumUnsortedVerticesC   -> vertices;   **int n = minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB
// compute the convex hull for minkowskiSumHull C       // **     minkowskiSumConvexHullListSizeC change = number of convex hull points
bool ContainsVertexA(double x, double y, int n) {
    for (int i = 0; i < n; i++) {
        if (minkowskiSumUnsortedVerticesA[i].x == x && minkowskiSumUnsortedVerticesA[i].y == y) return true;
    }
    return false;
}
void GeneratingVertexA(double xmin, double xmax, double ymin, double ymax) {
    for (int i = 0; i < 6; i++) { // x = [400, 1300], y = [100, 700]
        double x = GetRandomDouble(xmin, xmax);
        double y = GetRandomDouble(ymin, ymax);
        while (ContainsVertexA(x, y, i)) { // generate and no repeat vertex
            x = GetRandomDouble(xmin, xmax);
            y = GetRandomDouble(ymin, ymax);
        }
        minkowskiSumUnsortedVerticesA[i].x = x;
        minkowskiSumUnsortedVerticesA[i].y = y;
    }
}
bool ContainsVertexB(double x, double y, int n) {
    for (int i = 0; i < n; i++) {
        if (minkowskiSumUnsortedVerticesB[i].x == x && minkowskiSumUnsortedVerticesB[i].y == y) return true;
    }
    return false;
}
void GeneratingVertexB(double xmin, double xmax, double ymin, double ymax) {
    for (int i = 0; i < 6; i++) { // x = [400, 1300], y = [100, 700]
        double x = GetRandomDouble(xmin, xmax);
        double y = GetRandomDouble(ymin, ymax);
        while (ContainsVertexB(x, y, i)) { // generate and no repeat vertex
            x = GetRandomDouble(xmin, xmax);
            y = GetRandomDouble(ymin, ymax);
        }
        minkowskiSumUnsortedVerticesB[i].x = x;
        minkowskiSumUnsortedVerticesB[i].y = y;
    }
}

void MinkowskiSumGrahamScanCodeForHullC(POINT vertices[], int n) {
    for (int i = 0; i < n; i++) {
        minkowskiSumConvexHullPointsC[i].x = vertices[i].x;
        minkowskiSumConvexHullPointsC[i].y = vertices[i].y;
    }
    double minY = minkowskiSumConvexHullPointsC[0].y;
    int minIndex = 0;
    // find bottommost left vertex
    for (int i = 1; i < n; i++) {
        double currentY = minkowskiSumConvexHullPointsC[i].y;
        if ((currentY < minY) || (minY == currentY && minkowskiSumConvexHullPointsC[i].x < minkowskiSumConvexHullPointsC[minIndex].x)) {
            minY = currentY;
            minIndex = i;
        }
    }
    // the first vertex is the bottommost left vertex
    Swap(minkowskiSumConvexHullPointsC[0], minkowskiSumConvexHullPointsC[minIndex]);
    v0 = minkowskiSumConvexHullPointsC[0];                                        //     <- need a new global variable  POINT type   was named minkowskiSumHullCv0
    qsort(&minkowskiSumConvexHullPointsC[1], n - 1, sizeof(POINT), Compare);

    int j = 1;
    for (int i = 1; i < n; i++) {
        while (i < n - 1 && Direction(v0, minkowskiSumConvexHullPointsC[i], minkowskiSumConvexHullPointsC[i + 1]) == 0) {   // <- change minkowskiSumHullCv0
            i++;
        }
        minkowskiSumConvexHullPointsC[j] = minkowskiSumConvexHullPointsC[i];
        j++;
    }
    if (j < 3) return;
    stack<POINT> vertexStack;
    vertexStack.push(minkowskiSumConvexHullPointsC[0]);
    vertexStack.push(minkowskiSumConvexHullPointsC[1]);
    vertexStack.push(minkowskiSumConvexHullPointsC[2]);
    for (int i = 3; i < j; i++) {
        while (vertexStack.size() > 1 && Direction(SecondVertexInStack(vertexStack), vertexStack.top(), minkowskiSumConvexHullPointsC[i]) != -1) {
            vertexStack.pop();
        }
        vertexStack.push(minkowskiSumConvexHullPointsC[i]);
    }
    // empty the convexHullList in global if it contain elements
    while (minkowskiSumConvexHullListSizeC != 0) {                              //      <- need a new global variable  int type was named  minkowskiSumConvexHullListSizeC
        minkowskiSumConvexHullListSizeC--;                                      //      < change minkowskiSumConvexHullListSizeC
    }
    int counter = 0;
    while (!vertexStack.empty()) {
        POINT v = vertexStack.top();
        minkowskiSumConvexHullPointsC[counter] = v;                   //      <- need a new global variable  POINT[] type was name minkowskiSumConvexHullPointsC
        counter++;
        minkowskiSumConvexHullListSizeC++;                                      //      <- change minkowskiSumConvexHullListSizeC
        vertexStack.pop();
    }
    minkowskiSumConvexHullPointsC[counter] = minkowskiSumConvexHullPointsC[0];     //      <- change minkowskiSumConvexHullPointsC
}

// compute and store as unsortB by using SortedA and SortedB
// compute the x,y of A+B   // center of origin using min x,y and max x,y.
void MinkowskiSum(double minX, double maxX, double minY, double maxY) {
    // center of origin in the screen
    double centerX = (maxX - minX) / 2 + minX;
    double centerY = (maxY - minY) / 2 + minY;
    // find the unsorted a + b for a in A and b in B

    int cCounter = 0;
    for (int i = 0; i < minkowskiSumConvexHullListSizeA; i++) { // i for A
        for (int j = 0; j < minkowskiSumConvexHullListSizeB; j++) { // j for B
            double ax = (double)minkowskiSumConvexHullPointsA[i].x - centerX;
            double ay = (double)minkowskiSumConvexHullPointsA[i].y - centerY;
            double bx = (double)minkowskiSumConvexHullPointsB[j].x - centerX;
            double by = (double)minkowskiSumConvexHullPointsB[j].y - centerY;
            // computing a + b;
            double cx = ax + bx + centerX;
            double cy = ay + by + centerY;
            minkowskiSumUnsortedVerticesC[cCounter].x = (long)cx;
            minkowskiSumUnsortedVerticesC[cCounter].y = (long)cy;
            cCounter++;  // point to next array element;
        }
    }
    // 
}

// compute and store as unsortB by using SortedA and SortedB
// compute the x,y of A-B   // center of origin using min x,y and max x,y.
void MinkowskiDifference(double minX, double maxX, double minY, double maxY) {
    // center of origin in the screen
    double centerX = (maxX - minX) / 2 + minX;
    double centerY = (maxY - minY) / 2 + minY;
    // find the unsorted a + b for a in A and b in B

    int cCounter = 0;
    for (int i = 0; i < minkowskiSumConvexHullListSizeA; i++) { // i for A
        for (int j = 0; j < minkowskiSumConvexHullListSizeB; j++) { // j for B
            double ax = (double)minkowskiSumConvexHullPointsA[i].x - centerX;
            double ay = (double)minkowskiSumConvexHullPointsA[i].y - centerY;
            double bx = (double)minkowskiSumConvexHullPointsB[j].x - centerX;
            double by = (double)minkowskiSumConvexHullPointsB[j].y - centerY;
            // computing a + b;
            double cx = ax - bx + centerX;
            double cy = ay - by + centerY;
            minkowskiSumUnsortedVerticesC[cCounter].x = (long)cx;
            minkowskiSumUnsortedVerticesC[cCounter].y = (long)cy;
            cCounter++;  // point to next array element;
        }
    }
    // 
}

/*

    -------------------------------------------------- Minkowski Sum Algorithm's function ended --------------------------------------------------------------

*/



//    ----------------------------------------------------- GJK Algorithm's function ----------------------------------------------------------------



// true means 2 segment colliding.
bool IsGJKcolliding() {
    for (int i = 0; i < minkowskiSumConvexHullListSizeA; i++) {
        if (IsVertexInsideConvexHull(minkowskiSumConvexHullPointsB, minkowskiSumConvexHullPointsA[i], minkowskiSumConvexHullListSizeB)) return true;
    }
    for (int i = 0; i < minkowskiSumConvexHullListSizeB; i++) {
        if (IsVertexInsideConvexHull(minkowskiSumConvexHullPointsA, minkowskiSumConvexHullPointsB[i], minkowskiSumConvexHullListSizeA)) return true;
    }
    return false;
}



//    -------------------------------------------------- GJK Algorithm's function ended --------------------------------------------------------------



LRESULT CALLBACK WindowProcedure(HWND, UINT, WPARAM, LPARAM);
void AddControls(HWND);
void RengerBackground(HWND);
void DrawPoly(HDC, HPEN, int, POINT[]);
void DrawSinglePoint(HDC, HPEN, POINT);
void DrawPoint(HDC, HPEN, int, POINT[]);
void DrawLines(HDC);
HWND leftSide;

int WINAPI wWinMain(_In_ HINSTANCE hInstance,
    _In_opt_ HINSTANCE hPrevInstance,
    _In_ LPWSTR    lpCmdLine,
    _In_ int       nCmdShow)
{
    WNDCLASSW wc = { 0 };

    wc.hbrBackground = (HBRUSH)COLOR_WINDOW;
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hInstance = hInstance;
    wc.lpszClassName = L"myWindowClass";
    wc.lpfnWndProc = WindowProcedure;

    HBRUSH brush = CreateSolidBrush(RGB(0, 0, 0));
    wc.hbrBackground = brush;

    if (!RegisterClassW(&wc))
        return -1;

    CreateWindowW(L"myWindowClass", L"ConvexHullAlgorithms", WS_OVERLAPPEDWINDOW | WS_VISIBLE | WS_CLIPCHILDREN, 0, 0, 1400, 800, NULL, NULL, NULL, NULL);

    MSG msg = { 0 };

    while (GetMessage((&msg), NULL, NULL, NULL)) {

        TranslateMessage(&msg);
        DispatchMessage(&msg);

    }

    return 0;
}

LRESULT CALLBACK WindowProcedure(HWND hWnd, UINT msg, WPARAM wp, LPARAM lp) {

    HDC hdc;

    static POINT pt;

    static BOOL convexHULL = FALSE;
    static BOOL quickHULL = FALSE;
    static BOOL minkowskiSum = FALSE;
    static BOOL minkowskiDifference = FALSE;
    static BOOL caseGJK = FALSE;

    static BOOL DragPoly_convexHull = FALSE;
    static BOOL DragPoint_convexHull = FALSE;

    static BOOL DragPoint_quickHull = FALSE;
    static int index_quickHull = 0;

    static BOOL DragPoint_minkowskiSumA = FALSE;
    static BOOL DragPoint_minkowskiSumB = FALSE;
    static BOOL DragPoly_minkowskiSumA = FALSE;
    static BOOL DragPoly_minkowskiSumB = FALSE;
    static int index_minkowskiSum = 0;

    static BOOL DragPoint_minkowskiDifferenceA = FALSE;
    static BOOL DragPoint_minkowskiDifferenceB = FALSE;
    static BOOL DragPoly_minkowskiDifferenceA = FALSE;
    static BOOL DragPoly_minkowskiDifferenceB = FALSE;
    static int index_minkowskiDifference = 0;

    static HPEN pen;
    pen = CreatePen(PS_SOLID, 3, RGB(255, 255, 255));
    static HPEN circlePen;

    switch (msg)
    {
    case WM_COMMAND:
        switch (wp)
        {
        case POINT_CONVEX_HULL:
            hdc = GetDC(hWnd);

            RengerBackground(hWnd);

            srand(time(NULL));
            GeneratingVertex(400, 1300, 100, 700);
            GrahamScanCode(randomVertices, numberOfVertex);

            centerC.x = 850;
            centerC.y = 400;
            DrawPoly(hdc, pen, listSize + 1, convexHullPoints);
            DrawSinglePoint(hdc, CreatePen(PS_SOLID, 10, RGB(255, 0, 0)), centerC);

            ReleaseDC(hWnd, hdc);

            convexHULL = TRUE;
            quickHULL = FALSE;
            minkowskiSum = FALSE;
            minkowskiDifference = FALSE;
            caseGJK = FALSE;
            break;
        case QUICKHULL:
            hdc = GetDC(hWnd);

            RengerBackground(hWnd);

            srand(time(NULL));
            GeneratingVertex(400, 1300, 100, 700);
            GrahamScanCode(randomVertices, numberOfVertex);

            DrawPoly(hdc, pen, listSize + 1, convexHullPoints);
            DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 12, randomVertices);

            ReleaseDC(hWnd, hdc);

            quickHULL = TRUE;
            convexHULL = FALSE;
            minkowskiSum = FALSE;
            minkowskiDifference = FALSE;
            caseGJK = FALSE;
            break;
        case MINKOWSKI_DIFFERENCE:
            hdc = GetDC(hWnd);

            RengerBackground(hWnd);
            DrawLines(hdc);

            srand(time(NULL));
            GeneratingVertexA(850, 1000, 300, 350);
            MinkowskiSumGrahamScanCodeForHullA(minkowskiSumUnsortedVerticesA, 6);

            DrawPoly(hdc, pen, minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
            DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);



            GeneratingVertexB(1000, 1150, 250, 300);
            MinkowskiSumGrahamScanCodeForHullB(minkowskiSumUnsortedVerticesB, 6);

            DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
            DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);


            MinkowskiDifference(230, 1400, 0, 800);
            MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);

            DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);

            ReleaseDC(hWnd, hdc);

            minkowskiDifference = TRUE;
            minkowskiSum = FALSE;
            quickHULL = FALSE;
            convexHULL = FALSE;
            caseGJK = FALSE;
            break;
        case MINKOWSKI_SUM:
            hdc = GetDC(hWnd);

            RengerBackground(hWnd);
            DrawLines(hdc);

            srand(time(NULL));
            GeneratingVertexA(850, 1000, 300, 350);
            MinkowskiSumGrahamScanCodeForHullA(minkowskiSumUnsortedVerticesA, 6);

            DrawPoly(hdc, pen, minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
            DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);



            GeneratingVertexB(1000, 1150, 250, 300);
            MinkowskiSumGrahamScanCodeForHullB(minkowskiSumUnsortedVerticesB, 6);

            DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
            DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);



            MinkowskiSum(230, 1400, 0, 800);
            MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);

            DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);

            ReleaseDC(hWnd, hdc);

            minkowskiSum = TRUE;
            quickHULL = FALSE;
            convexHULL = FALSE;
            minkowskiDifference = FALSE;
            caseGJK = FALSE;
            break;
        case GJK:
            hdc = GetDC(hWnd);

            RengerBackground(hWnd);
            DrawLines(hdc);

            srand(time(NULL));
            GeneratingVertexA(850, 1000, 300, 350);
            MinkowskiSumGrahamScanCodeForHullA(minkowskiSumUnsortedVerticesA, 6);

            DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 0)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
            DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);



            GeneratingVertexB(1000, 1150, 250, 300);
            MinkowskiSumGrahamScanCodeForHullB(minkowskiSumUnsortedVerticesB, 6);

            DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(0, 0, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
            DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);

            MinkowskiDifference(230, 1400, 0, 800);
            MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);

            if (IsGJKcolliding()) {
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(0, 128, 0)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
            }
            else {
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
            }

            ReleaseDC(hWnd, hdc);

            caseGJK = TRUE;
            minkowskiDifference = FALSE;
            minkowskiSum = FALSE;
            quickHULL = FALSE;
            convexHULL = FALSE;
            break;
        }
        break;
    case WM_CREATE:
        AddControls(hWnd);
        break;
    case WM_LBUTTONDOWN:
        pt.x = (LONG)LOWORD(lp);
        pt.y = (LONG)HIWORD(lp);
        if (convexHULL) {
            if (((pt.x - centerC.x) <= 5 && (pt.x - centerC.x) >= -5) && ((pt.y - centerC.y) <= 5 && (pt.y - centerC.y) >= -5)) {
                DragPoint_convexHull = TRUE;
            }
            else if (IsVertexInsideConvexHull(convexHullPoints, pt, listSize)) {
                DragPoly_convexHull = TRUE;
            }
        }
        if (quickHULL) {
            int counter = 0;
            for (int i = 0; i < 12; i++) {
                if (((pt.x - randomVertices[i].x) <= 5 && (pt.x - randomVertices[i].x) >= -5) && ((pt.y - randomVertices[i].y) <= 5 && (pt.y - randomVertices[i].y) >= -5)) {
                    index_quickHull = counter;
                    DragPoint_quickHull = TRUE;
                    break;
                }
                counter++;
            }
        }
        if (minkowskiSum) {

            int counter = 0;
            for (int i = 0; i < 6; i++) {
                if (((pt.x - minkowskiSumUnsortedVerticesA[i].x) <= 3 && (pt.x - minkowskiSumUnsortedVerticesA[i].x) >= -3) && ((pt.y - minkowskiSumUnsortedVerticesA[i].y) <= 3 && (pt.y - minkowskiSumUnsortedVerticesA[i].y) >= -3)) {
                    index_minkowskiSum = counter;
                    DragPoint_minkowskiSumA = TRUE;
                    break;
                }
                counter++;
            }
            if (!DragPoint_minkowskiSumA) {
                counter = 0;
                for (int i = 0; i < 6; i++) {
                    if (((pt.x - minkowskiSumUnsortedVerticesB[i].x) <= 3 && (pt.x - minkowskiSumUnsortedVerticesB[i].x) >= -3) && ((pt.y - minkowskiSumUnsortedVerticesB[i].y) <= 3 && (pt.y - minkowskiSumUnsortedVerticesB[i].y) >= -3)) {
                        index_minkowskiSum = counter;
                        DragPoint_minkowskiSumB = TRUE;
                        break;
                    }
                    counter++;
                }
            }
            if (!DragPoint_minkowskiSumA && !DragPoint_minkowskiSumB && IsVertexInsideConvexHull(minkowskiSumConvexHullPointsA, pt, minkowskiSumConvexHullListSizeA) == 1) {
                DragPoly_minkowskiSumA = TRUE;
            }
            if (!DragPoint_minkowskiSumA && !DragPoint_minkowskiSumB && IsVertexInsideConvexHull(minkowskiSumConvexHullPointsB, pt, minkowskiSumConvexHullListSizeB) == 1) {
                DragPoly_minkowskiSumB = TRUE;
            }
        }
        if (minkowskiDifference) {

            int counter = 0;
            for (int i = 0; i < 6; i++) {
                if (((pt.x - minkowskiSumUnsortedVerticesA[i].x) <= 3 && (pt.x - minkowskiSumUnsortedVerticesA[i].x) >= -3) && ((pt.y - minkowskiSumUnsortedVerticesA[i].y) <= 3 && (pt.y - minkowskiSumUnsortedVerticesA[i].y) >= -3)) {
                    index_minkowskiDifference = counter;
                    DragPoint_minkowskiDifferenceA = TRUE;
                    break;
                }
                counter++;
            }
            if (!DragPoint_minkowskiDifferenceA) {
                counter = 0;
                for (int i = 0; i < 6; i++) {
                    if (((pt.x - minkowskiSumUnsortedVerticesB[i].x) <= 3 && (pt.x - minkowskiSumUnsortedVerticesB[i].x) >= -3) && ((pt.y - minkowskiSumUnsortedVerticesB[i].y) <= 3 && (pt.y - minkowskiSumUnsortedVerticesB[i].y) >= -3)) {
                        index_minkowskiDifference = counter;
                        DragPoint_minkowskiDifferenceB = TRUE;
                        break;
                    }
                    counter++;
                }
            }
            if (!DragPoint_minkowskiDifferenceA && !DragPoint_minkowskiDifferenceB && IsVertexInsideConvexHull(minkowskiSumConvexHullPointsA, pt, minkowskiSumConvexHullListSizeA) == 1) {
                DragPoly_minkowskiDifferenceA = TRUE;
            }
            if (!DragPoint_minkowskiDifferenceA && !DragPoint_minkowskiDifferenceB && IsVertexInsideConvexHull(minkowskiSumConvexHullPointsB, pt, minkowskiSumConvexHullListSizeB) == 1) {
                DragPoly_minkowskiDifferenceB = TRUE;
            }
        }


        if (caseGJK) {

            int counter = 0;
            for (int i = 0; i < 6; i++) {
                if (((pt.x - minkowskiSumUnsortedVerticesA[i].x) <= 3 && (pt.x - minkowskiSumUnsortedVerticesA[i].x) >= -3) && ((pt.y - minkowskiSumUnsortedVerticesA[i].y) <= 3 && (pt.y - minkowskiSumUnsortedVerticesA[i].y) >= -3)) {
                    index_minkowskiDifference = counter;
                    DragPoint_minkowskiDifferenceA = TRUE;
                    break;
                }
                counter++;
            }
            if (!DragPoint_minkowskiDifferenceA) {
                counter = 0;
                for (int i = 0; i < 6; i++) {
                    if (((pt.x - minkowskiSumUnsortedVerticesB[i].x) <= 3 && (pt.x - minkowskiSumUnsortedVerticesB[i].x) >= -3) && ((pt.y - minkowskiSumUnsortedVerticesB[i].y) <= 3 && (pt.y - minkowskiSumUnsortedVerticesB[i].y) >= -3)) {
                        index_minkowskiDifference = counter;
                        DragPoint_minkowskiDifferenceB = TRUE;
                        break;
                    }
                    counter++;
                }
            }
            if (!DragPoint_minkowskiDifferenceA && !DragPoint_minkowskiDifferenceB && IsVertexInsideConvexHull(minkowskiSumConvexHullPointsA, pt, minkowskiSumConvexHullListSizeA) == 1) {
                DragPoly_minkowskiDifferenceA = TRUE;
            }
            if (!DragPoint_minkowskiDifferenceA && !DragPoint_minkowskiDifferenceB && IsVertexInsideConvexHull(minkowskiSumConvexHullPointsB, pt, minkowskiSumConvexHullListSizeB) == 1) {
                DragPoly_minkowskiDifferenceB = TRUE;
            }
        }

        break;
    case WM_MOUSEMOVE:
        hdc = GetDC(hWnd);

        if (convexHULL) {
            if (DragPoly_convexHull) {
                RengerBackground(hWnd);

                for (int i = 0; i < listSize + 1; i++) {
                    convexHullPoints[i].x += (LOWORD(lp) - pt.x);
                    convexHullPoints[i].y += (HIWORD(lp) - pt.y);
                }

                DrawPoly(hdc, pen, listSize + 1, convexHullPoints);

                if (IsVertexInsideConvexHull(convexHullPoints, centerC, listSize)) {
                    circlePen = CreatePen(PS_SOLID, 10, RGB(255, 0, 0));
                }
                else {
                    circlePen = CreatePen(PS_SOLID, 10, RGB(50, 205, 50));
                }
                DrawSinglePoint(hdc, circlePen, centerC);
            }

            if (DragPoint_convexHull) {
                RengerBackground(hWnd);

                if (IsVertexInsideConvexHull(convexHullPoints, centerC, listSize)) {
                    circlePen = CreatePen(PS_SOLID, 10, RGB(255, 0, 0));
                }
                else {
                    circlePen = CreatePen(PS_SOLID, 10, RGB(50, 205, 50));
                }

                centerC.x += (LOWORD(lp) - centerC.x);
                centerC.y += (HIWORD(lp) - centerC.y);
                DrawSinglePoint(hdc, circlePen, centerC);

                DrawPoly(hdc, pen, listSize + 1, convexHullPoints);
            }
        }


        if (quickHULL) {
            if (DragPoint_quickHull) {
                RengerBackground(hWnd);

                randomVertices[index_quickHull].x += (LOWORD(lp) - randomVertices[index_quickHull].x);
                randomVertices[index_quickHull].y += (HIWORD(lp) - randomVertices[index_quickHull].y);

                GrahamScanCode(randomVertices, 12);

                DrawPoly(hdc, pen, listSize + 1, convexHullPoints);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 12, randomVertices);
            }
        }


        if (minkowskiSum) {
            if (DragPoint_minkowskiSumA) {
                RengerBackground(hWnd);
                DrawLines(hdc);

                minkowskiSumUnsortedVerticesA[index_minkowskiSum].x += (LOWORD(lp) - minkowskiSumUnsortedVerticesA[index_minkowskiSum].x);
                minkowskiSumUnsortedVerticesA[index_minkowskiSum].y += (HIWORD(lp) - minkowskiSumUnsortedVerticesA[index_minkowskiSum].y);

                MinkowskiSumGrahamScanCodeForHullA(minkowskiSumUnsortedVerticesA, 6);

                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);

                MinkowskiSum(230, 1400, 0, 800);
                MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
            }
            else if (DragPoint_minkowskiSumB) {
                RengerBackground(hWnd);
                DrawLines(hdc);

                minkowskiSumUnsortedVerticesB[index_minkowskiSum].x += (LOWORD(lp) - minkowskiSumUnsortedVerticesB[index_minkowskiSum].x);
                minkowskiSumUnsortedVerticesB[index_minkowskiSum].y += (HIWORD(lp) - minkowskiSumUnsortedVerticesB[index_minkowskiSum].y);

                MinkowskiSumGrahamScanCodeForHullB(minkowskiSumUnsortedVerticesB, 6);

                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);

                MinkowskiSum(230, 1400, 0, 800);
                MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
            }
            else if (DragPoly_minkowskiSumA) {
                RengerBackground(hWnd);
                DrawLines(hdc);
                for (int i = 0; i < 6; i++) {
                    minkowskiSumUnsortedVerticesA[i].x += (LOWORD(lp) - pt.x);
                    minkowskiSumUnsortedVerticesA[i].y += (HIWORD(lp) - pt.y);
                }
                MinkowskiSumGrahamScanCodeForHullA(minkowskiSumUnsortedVerticesA, 6);

                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);

                MinkowskiSum(230, 1400, 0, 800);
                MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);

            }
            else if (DragPoly_minkowskiSumB) {
                RengerBackground(hWnd);
                DrawLines(hdc);
                for (int i = 0; i < 6; i++) {
                    minkowskiSumUnsortedVerticesB[i].x += (LOWORD(lp) - pt.x);
                    minkowskiSumUnsortedVerticesB[i].y += (HIWORD(lp) - pt.y);
                }
                MinkowskiSumGrahamScanCodeForHullB(minkowskiSumUnsortedVerticesB, 6);

                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);

                MinkowskiSum(230, 1400, 0, 800);
                MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);

            }
        }

        if (minkowskiDifference) {
            if (DragPoint_minkowskiDifferenceA) {
                RengerBackground(hWnd);
                DrawLines(hdc);

                minkowskiSumUnsortedVerticesA[index_minkowskiDifference].x += (LOWORD(lp) - minkowskiSumUnsortedVerticesA[index_minkowskiDifference].x);
                minkowskiSumUnsortedVerticesA[index_minkowskiDifference].y += (HIWORD(lp) - minkowskiSumUnsortedVerticesA[index_minkowskiDifference].y);

                MinkowskiSumGrahamScanCodeForHullA(minkowskiSumUnsortedVerticesA, 6);

                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);

                MinkowskiDifference(230, 1400, 0, 800);
                MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
            }
            else if (DragPoint_minkowskiDifferenceB) {
                RengerBackground(hWnd);
                DrawLines(hdc);

                minkowskiSumUnsortedVerticesB[index_minkowskiDifference].x += (LOWORD(lp) - minkowskiSumUnsortedVerticesB[index_minkowskiDifference].x);
                minkowskiSumUnsortedVerticesB[index_minkowskiDifference].y += (HIWORD(lp) - minkowskiSumUnsortedVerticesB[index_minkowskiDifference].y);

                MinkowskiSumGrahamScanCodeForHullB(minkowskiSumUnsortedVerticesB, 6);

                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);

                MinkowskiDifference(230, 1400, 0, 800);
                MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
            }
            else if (DragPoly_minkowskiDifferenceA) {
                RengerBackground(hWnd);
                DrawLines(hdc);
                for (int i = 0; i < 6; i++) {
                    minkowskiSumUnsortedVerticesA[i].x += (LOWORD(lp) - pt.x);
                    minkowskiSumUnsortedVerticesA[i].y += (HIWORD(lp) - pt.y);
                }
                MinkowskiSumGrahamScanCodeForHullA(minkowskiSumUnsortedVerticesA, 6);

                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);

                MinkowskiDifference(230, 1400, 0, 800);
                MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);

            }
            else if (DragPoly_minkowskiDifferenceB) {
                RengerBackground(hWnd);
                DrawLines(hdc);
                for (int i = 0; i < 6; i++) {
                    minkowskiSumUnsortedVerticesB[i].x += (LOWORD(lp) - pt.x);
                    minkowskiSumUnsortedVerticesB[i].y += (HIWORD(lp) - pt.y);
                }
                MinkowskiSumGrahamScanCodeForHullB(minkowskiSumUnsortedVerticesB, 6);

                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);

                MinkowskiDifference(230, 1400, 0, 800);
                MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);

            }
        }


        if (caseGJK) {
            if (DragPoint_minkowskiDifferenceA) {
                RengerBackground(hWnd);
                DrawLines(hdc);

                minkowskiSumUnsortedVerticesA[index_minkowskiDifference].x += (LOWORD(lp) - minkowskiSumUnsortedVerticesA[index_minkowskiDifference].x);
                minkowskiSumUnsortedVerticesA[index_minkowskiDifference].y += (HIWORD(lp) - minkowskiSumUnsortedVerticesA[index_minkowskiDifference].y);

                MinkowskiSumGrahamScanCodeForHullA(minkowskiSumUnsortedVerticesA, 6);

                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 0)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(0, 0, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);

                MinkowskiDifference(230, 1400, 0, 800);
                MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);
                if (IsGJKcolliding()) {
                    DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(0, 128, 0)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                }
                else {
                    DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                }
            }
            else if (DragPoint_minkowskiDifferenceB) {
                RengerBackground(hWnd);
                DrawLines(hdc);

                minkowskiSumUnsortedVerticesB[index_minkowskiDifference].x += (LOWORD(lp) - minkowskiSumUnsortedVerticesB[index_minkowskiDifference].x);
                minkowskiSumUnsortedVerticesB[index_minkowskiDifference].y += (HIWORD(lp) - minkowskiSumUnsortedVerticesB[index_minkowskiDifference].y);

                MinkowskiSumGrahamScanCodeForHullB(minkowskiSumUnsortedVerticesB, 6);

                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 0)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(0, 0, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);

                MinkowskiDifference(230, 1400, 0, 800);
                MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);
                if (IsGJKcolliding()) {
                    DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(0, 128, 0)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                }
                else {
                    DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                }
            }
            else if (DragPoly_minkowskiDifferenceA) {
                RengerBackground(hWnd);
                DrawLines(hdc);
                for (int i = 0; i < 6; i++) {
                    minkowskiSumUnsortedVerticesA[i].x += (LOWORD(lp) - pt.x);
                    minkowskiSumUnsortedVerticesA[i].y += (HIWORD(lp) - pt.y);
                }
                MinkowskiSumGrahamScanCodeForHullA(minkowskiSumUnsortedVerticesA, 6);

                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 0)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(0, 0, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);

                MinkowskiDifference(230, 1400, 0, 800);
                MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);
                if (IsGJKcolliding()) {
                    DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(0, 128, 0)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                }
                else {
                    DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                }

            }
            else if (DragPoly_minkowskiDifferenceB) {
                RengerBackground(hWnd);
                DrawLines(hdc);
                for (int i = 0; i < 6; i++) {
                    minkowskiSumUnsortedVerticesB[i].x += (LOWORD(lp) - pt.x);
                    minkowskiSumUnsortedVerticesB[i].y += (HIWORD(lp) - pt.y);
                }
                MinkowskiSumGrahamScanCodeForHullB(minkowskiSumUnsortedVerticesB, 6);

                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 0)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(0, 0, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), 6, minkowskiSumUnsortedVerticesB);

                MinkowskiDifference(230, 1400, 0, 800);
                MinkowskiSumGrahamScanCodeForHullC(minkowskiSumUnsortedVerticesC, minkowskiSumConvexHullListSizeA * minkowskiSumConvexHullListSizeB);
                if (IsGJKcolliding()) {
                    DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(0, 128, 0)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                }
                else {
                    DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                }

            }
        }


        pt.x = (LONG)LOWORD(lp);
        pt.y = (LONG)HIWORD(lp);
        ReleaseDC(hWnd, hdc);
        break;
    case WM_LBUTTONUP:
        hdc = GetDC(hWnd);


        if (convexHULL) {
            if (DragPoly_convexHull)
            {
                DrawPoly(hdc, pen, listSize + 1, convexHullPoints);
                DragPoly_convexHull = FALSE;
            }
            if (DragPoint_convexHull) {
                DrawSinglePoint(hdc, circlePen, centerC);
                DragPoint_convexHull = FALSE;
            }
        }


        if (quickHULL) {
            if (DragPoint_quickHull) {
                DrawSinglePoint(hdc, circlePen, randomVertices[index_quickHull]);
                DragPoint_quickHull = FALSE;
            }
        }


        if (minkowskiSum) {
            if (DragPoint_minkowskiSumA) {
                DrawSinglePoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), minkowskiSumUnsortedVerticesA[index_minkowskiSum]);
                DragPoint_minkowskiSumA = FALSE;
            }
            else if (DragPoint_minkowskiSumB) {
                DrawSinglePoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), minkowskiSumUnsortedVerticesB[index_minkowskiSum]);
                DragPoint_minkowskiSumB = FALSE;
            }
            else if (DragPoly_minkowskiSumA) {
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                DragPoly_minkowskiSumA = FALSE;
            }
            else if (DragPoly_minkowskiSumB) {
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                DragPoly_minkowskiSumB = FALSE;
            }
        }

        if (minkowskiDifference) {
            if (DragPoint_minkowskiDifferenceA) {
                DrawSinglePoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), minkowskiSumUnsortedVerticesA[index_minkowskiSum]);
                DragPoint_minkowskiDifferenceA = FALSE;
            }
            else if (DragPoint_minkowskiDifferenceB) {
                DrawSinglePoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), minkowskiSumUnsortedVerticesB[index_minkowskiSum]);
                DragPoint_minkowskiDifferenceB = FALSE;
            }
            else if (DragPoly_minkowskiDifferenceA) {
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                DragPoly_minkowskiDifferenceA = FALSE;
            }
            else if (DragPoly_minkowskiDifferenceB) {
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 255, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                DragPoly_minkowskiDifferenceB = FALSE;
            }
        }


        if (caseGJK) {
            if (DragPoint_minkowskiDifferenceA) {
                DrawSinglePoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), minkowskiSumUnsortedVerticesA[index_minkowskiSum]);
                DragPoint_minkowskiDifferenceA = FALSE;
            }
            else if (DragPoint_minkowskiDifferenceB) {
                DrawSinglePoint(hdc, CreatePen(PS_SOLID, 10, RGB(50, 205, 50)), minkowskiSumUnsortedVerticesB[index_minkowskiSum]);
                DragPoint_minkowskiDifferenceB = FALSE;
            }
            else if (DragPoly_minkowskiDifferenceA) {
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 0)), minkowskiSumConvexHullListSizeA + 1, minkowskiSumConvexHullPointsA);
                if (IsGJKcolliding()) {
                    DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(0, 128, 0)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                }
                else {
                    DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                }
                DragPoly_minkowskiDifferenceA = FALSE;
            }
            else if (DragPoly_minkowskiDifferenceB) {
                DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(0, 0, 255)), minkowskiSumConvexHullListSizeB + 1, minkowskiSumConvexHullPointsB);
                if (IsGJKcolliding()) {
                    DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(0, 128, 0)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                }
                else {
                    DrawPoly(hdc, CreatePen(PS_SOLID, 3, RGB(255, 0, 255)), minkowskiSumConvexHullListSizeC + 1, minkowskiSumConvexHullPointsC);
                }
                DragPoly_minkowskiDifferenceB = FALSE;
            }
        }

        ReleaseDC(hWnd, hdc);
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    case WM_CTLCOLORSTATIC:
        return (INT_PTR)CreateSolidBrush(RGB(0, 0, 0));
    default:
        return DefWindowProcW(hWnd, msg, wp, lp);
    }

}


void AddControls(HWND hWnd)
{
    HWND leftSide = CreateWindow(L"Static", L"", WS_VISIBLE | WS_CHILD, 0, 0, 230, 850, hWnd, NULL, NULL, NULL);
    CreateWindowW(L"Button", L"Point Convex Hull", WS_VISIBLE | WS_CHILD, 20, 20, 200, 40, hWnd, (HMENU)POINT_CONVEX_HULL, NULL, NULL);
    CreateWindowW(L"Button", L"Quickhull", WS_VISIBLE | WS_CHILD, 20, 80, 200, 40, hWnd, (HMENU)QUICKHULL, NULL, NULL);
    CreateWindowW(L"Button", L"Minkowski Difference", WS_VISIBLE | WS_CHILD, 20, 140, 200, 40, hWnd, (HMENU)MINKOWSKI_DIFFERENCE, NULL, NULL);
    CreateWindowW(L"Button", L"Minkowski Sum", WS_VISIBLE | WS_CHILD, 20, 200, 200, 40, hWnd, (HMENU)MINKOWSKI_SUM, NULL, NULL);
    CreateWindowW(L"Button", L"GJK", WS_VISIBLE | WS_CHILD, 20, 260, 200, 40, hWnd, (HMENU)GJK, NULL, NULL);
}

void RengerBackground(HWND hWnd) {
    HDC hdc = GetDC(hWnd);
    HPEN hpen_background = CreatePen(PS_SOLID, 7, RGB(255, 255, 255));
    HBRUSH hbrush_background = CreateSolidBrush(RGB(0, 0, 0));

    SelectObject(hdc, hpen_background);
    SelectObject(hdc, hbrush_background);

    Rectangle(hdc, 230, 0, 1550, 850);

    DeleteObject(hpen_background);
    DeleteObject(hbrush_background);

    ReleaseDC(hWnd, hdc);
}

void DrawPoly(HDC hdc, HPEN pen, int size, POINT vertices[]) {
    SelectObject(hdc, pen);
    Polyline(hdc, vertices, size);
    DeleteObject(pen);
}

void DrawPoint(HDC hdc, HPEN pen, int size, POINT vertices[]) {
    SelectObject(hdc, pen);
    for (int i = 0; i < size; i++) {
        Ellipse(hdc, vertices[i].x, vertices[i].y, vertices[i].x + 5, vertices[i].y + 5);
    }
    DeleteObject(pen);
};

void DrawSinglePoint(HDC hdc, HPEN pen, POINT point) {
    SelectObject(hdc, pen);
    Ellipse(hdc, point.x, point.y, point.x + 5, point.y + 5);
    DeleteObject(pen);
};

void DrawLines(HDC hdc) {
    HPEN background = CreatePen(PS_SOLID, 1, RGB(128, 128, 128));
    SelectObject(hdc, background);
    int xAxis = 0;
    int yAxis = 0;
    for (int i = 0; i < 150; i++) {
        xAxis = i * 15;
        MoveToEx(hdc, xAxis, 0, NULL);
        LineTo(hdc, xAxis, 2000);
    }
    for (int i = 0; i < 100; i++) {
        yAxis = i * 15;
        MoveToEx(hdc, 0, yAxis, NULL);
        LineTo(hdc, 2000, yAxis);
    }
    DeleteObject(background);

    background = CreatePen(PS_SOLID, 3, RGB(128, 128, 128));
    SelectObject(hdc, background);
    MoveToEx(hdc, 838, 0, NULL);
    LineTo(hdc, 838, 2000);
    MoveToEx(hdc, 0, 390, NULL);
    LineTo(hdc, 2000, 390);
    DeleteObject(background);
};
