#include <bits/stdc++.h>
#include <iostream>

using namespace std;

#define PI acos(-1.0)

ofstream stage1_out;
ofstream stage2_out;
ofstream stage3_out;

double eyeX, eyeY, eyeZ;
double lookX, lookY, lookZ;
double upX, upY, upZ;
double fovY, aspectRatio, near, far;
double tx, ty, tz;
double sx, sy, sz;
double angle, ax, ay, az;
double p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z;

double lx, ly, lz; // for stage 2
double rx, ry, rz;
double ux, uy, uz;
double L;
double **V, **R, **T;

double fovX, t, r; //for stage 3

double **P;

stack<pair<double**, bool>> stack_stack;

double** zero_matrix()
{
    double** matrix = new double*[4];
    for(int i = 0; i < 4; i++)
    {
        matrix[i] = new double[4];
    }

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            matrix[i][j] = 0;
        }
    }
    return matrix;
}

double** multiply(double** matrix1, double** matrix2)
{
    double** matrix = zero_matrix();

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                matrix[i][j] += matrix1[i][k] * matrix2[k][j];

    return matrix;
}

double** I()
{
    double** identity_matrix = zero_matrix();
 
    for(int i = 0; i < 4; i++)
    {
        identity_matrix[i][i] = 1.0;
    }

    return identity_matrix;
}

//multiplies the given transformation matrix with stack_stack.top() and pushes the product back into stack_stack.
void multiply_transformation(double** matrix)
{
    double** prev = stack_stack.top().first;
    double** result = multiply(prev, matrix);

    stack_stack.push(make_pair(result, false));
}

// when pop() is called, is_pushed is checked
void pop()
{
    while(stack_stack.top().second == false)
    {
        stack_stack.pop();
    }
    stack_stack.pop();
}

// pushes a copy of stack_stack.top() with is_pushed set to true
void push()
{
    double** matrix = stack_stack.top().first;
    stack_stack.push(make_pair(matrix, true));
}

void rotate(double angle, double ax, double ay, double az)
{
    double cos_theta = cos(PI * angle / 180.0);
    double sin_theta = sin(PI * angle / 180.0);

    // normalize the axis of rotation
    double L = sqrt(ax * ax + ay * ay + az * az);
    double ux = ax / L;
    double uy = ay / L;
    double uz = az / L;

    double** matrix = zero_matrix();

    matrix[0][0] = cos_theta + ux * ux * (1 - cos_theta);
    matrix[0][1] = ux * uy * (1 - cos_theta) - uz * sin_theta;
    matrix[0][2] = ux * uz * (1 - cos_theta) + uy * sin_theta;

    matrix[1][0] = ux * uy * (1 - cos_theta) + uz * sin_theta;
    matrix[1][1] = cos_theta + uy * uy * (1 - cos_theta);
    matrix[1][2] = uy * uz * (1 - cos_theta) - ux * sin_theta;

    matrix[2][0] = ux * uz * (1 - cos_theta) - uy * sin_theta;
    matrix[2][1] = uy * uz * (1 - cos_theta) + ux * sin_theta;
    matrix[2][2] = cos_theta + uz * uz * (1 - cos_theta);

    matrix[3][3] = 1.0;

    multiply_transformation(matrix);
}

void scale(double sx, double sy, double sz)
{
    double** matrix = I();

    matrix[0][0] = sx;
    matrix[1][1] = sy;
    matrix[2][2] = sz;

    multiply_transformation(matrix);
}

void translate(double tx, double ty, double tz)
{
    double** matrix = I();

    matrix[0][3] = tx;
    matrix[1][3] = ty;
    matrix[2][3] = tz;

    multiply_transformation(matrix);
}

double** triangle(double x0, double y0, double z0, 
            double x1, double y1, double z1,
            double x2, double y2, double z2)
{
    double** matrix = zero_matrix();

    matrix[0][0] = x0;
    matrix[1][0] = y0;
    matrix[2][0] = z0;
    matrix[3][0] = 1.0;

    matrix[0][1] = x1;
    matrix[1][1] = y1;
    matrix[2][1] = z1;
    matrix[3][1] = 1.0;

    matrix[0][2] = x2;
    matrix[1][2] = y2;
    matrix[2][2] = z2;
    matrix[3][2] = 1.0;

    matrix[0][3] = matrix[1][3] = matrix[2][3] = matrix[3][3] = 1.0;

    double** result = multiply(stack_stack.top().first, matrix);
    return result;
}


void print_triangle(ofstream &stream, double** triangle)
{
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                stream << fixed << setprecision(7) << triangle[j][i] << " ";
            }
            stream << endl;
        }
        stream << endl;
}

void init()
{
    stack_stack.push(make_pair(I(), false));

    stage1_out.open("stage1.txt");
    stage2_out.open("stage2.txt");
    stage3_out.open("stage3.txt");

    cin >> eyeX >> eyeY >> eyeZ;
    cin >> lookX >> lookY >> lookZ;
    cin >> upX >> upY >> upZ;
    cin >> fovY >> aspectRatio >> near >> far;

    lx = lookX - eyeX;
    ly = lookY - eyeY;
    lz = lookZ - eyeZ;

    L = sqrt(lx * lx + ly * ly + lz * lz);

    lx = lx / L;
    ly = ly / L;
    lz = lz / L;

    rx = ly * upZ - lz * upY;
    ry = lz * upX - lx * upZ;
    rz = lx * upY - ly * upX;

    L = sqrt(rx * rx + ry * ry + rz * rz);

    rx = rx / L;
    ry = ry / L;
    rz = rz / L;

    ux = ry * lz - rz * ly;
    uy = rz * lx - rx * lz;
    uz = rx * ly - ry * lx;

    L = sqrt(ux * ux + uy * uy + uz * uz);

    ux = ux / L;
    uy = uy / L;
    uz = uz / L;

    R = zero_matrix();
    
    R[0][0] = rx;
    R[0][1] = ry;
    R[0][2] = rz;
    R[1][0] = ux;
    R[1][1] = uy;
    R[1][2] = uz;
    R[2][0] = -lx;
    R[2][1] = -ly;
    R[2][2] = -lz;
    R[3][3] = 1.0;

    T = I();

    T[0][3] = -eyeX;
    T[1][3] = -eyeY;
    T[2][3] = -eyeZ;

    V = multiply(R, T);

    fovX = fovY * aspectRatio; //stage 3
    t = near * tan(fovY /2.0 * PI / 180.0);
    r = near * tan(fovX /2.0 * PI / 180.0);

    P = zero_matrix();

    P[0][0] = near / r;
    P[1][1] = near / t;
    P[2][2] = -(far + near) / (far - near);
    P[2][3] = -(2 * far * near) / (far - near);
    P[3][2] = -1.0;
}

int main()
{
    string command;
    freopen("scene.txt", "r", stdin);
    init();

    while (1)
    {
        cin >> command;

        if (command == "triangle")
        {
            cin >> p1x >> p1y >> p1z;
            cin >> p2x >> p2y >> p2z;
            cin >> p3x >> p3y >> p3z;

            double **triangle1 = triangle(p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    if(j==2)
                    {
                        stage1_out << fixed << setprecision(7) << triangle1[j][i];
                    }
                    else
                    {
                        stage1_out << fixed << setprecision(7) << triangle1[j][i] << " ";
                    }
                }
                stage1_out << endl;
            }
            stage1_out << endl;

            double **triangle2 = multiply(V, triangle1);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    if(j==2)
                    {
                        stage2_out << fixed << setprecision(7) << triangle2[j][i];
                    }
                    else
                    {
                        stage2_out << fixed << setprecision(7) << triangle2[j][i] << " ";
                    }
                }
                stage2_out << endl;
            }
            stage2_out << endl;

            double **triangle3 = multiply(P, triangle2);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    if(j==2)
                    {
                        stage3_out << fixed << setprecision(7) << triangle3[j][i]/triangle3[3][i];
                    }
                    else
                    {
                        stage3_out << fixed << setprecision(7) << triangle3[j][i]/triangle3[3][i] << " ";
                    }
                    
                }
                stage3_out << endl;
            }
            stage3_out << endl;
        }
        else if (command == "translate")
        {
            cin >> tx >> ty >> tz;
            translate(tx, ty, tz);
        }
        else if (command == "scale")
        {
            cin >> sx >> sy >> sz;
            scale(sx, sy, sz);
        }
        else if (command == "rotate")
        {
            cin >> angle >> ax >> ay >> az;
            rotate(angle, ax, ay, az);
        }
        else if (command == "push")
        {
            push();
        }
        else if (command == "pop")
        {
            pop();
        }
        else if (command == "end")
        {
            break;
        }
        else
        {
            cout << "wrong command" << endl;
            break;
        }
    }
    return 0;
}