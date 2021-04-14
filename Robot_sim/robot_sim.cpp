#include <iostream>
#include <GL/glut.h>
#include <vector>
#include <stdarg.h>
#include <fstream>
#include <ctime>
#include "omp.h"
#include <algorithm>
#include <math.h>
#include <cmath>
#include <numeric>




#define LEN 2000  

#define Cos(th) cos(PI/180*(th))
#define Sin(th) sin(PI/180*(th))

#define PI 3.1412


struct Mass
{
    double m;      
    double pos[3];    
    double vel[3];    
    double acc[3];    
};

struct Spring
{
    double k;       
    double L_0;   
    int m1;        
    int m2;         
};

struct Gene
{
  
    double A;
    double B;
    double C;
};




int th = 0;            
int ph = 10;         

double asp = 1;    
int fov = 45;        
double dim = 10;  





int generationNumber = 2;
int cubotNumber = 10;
int simulationTime = 600;



unsigned int grassTexture;
unsigned int slimeTexture;
unsigned int skyBoxTexture[10];


double mass = 0.1;
double length = 1;
double gravity = 9.81;
double T = 0;

double timeStep = 0.001;
double restoreConstant = 100000;
double springConstant = 10000;
double dampingConstant = 0.998;
double frictionCoefficient = 0.8;

static GLint Frames = 0;
static GLfloat fps = -1;
static GLint T0 = 0;

GLfloat worldRotation[16] = { 1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1 };

std::ofstream bestGene;
std::ofstream popDis;


double calcMag(double x[], std::size_t sz)
{
    return std::sqrt(std::inner_product(x, x + sz, x, 0.0));
}

std::vector<int> sort(const std::vector<double>& v) {

    
    std::vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    
    sort(idx.begin(), idx.end(),
        [&v](int i1, int i2) {return v[i1] > v[i2]; });
    return idx;
}

class cubot
{
private:
    std::vector<Gene> gene;
public:
    double initialLocation[3] = { 0,0,0 };
    std::vector<Mass> cubemasses;
    std::vector<Spring> cubesprings;

    cubot(double initialX, double initialY, double initialZ, std::vector<Gene> Gene)
    {
        
        initialLocation[0] = initialX; initialLocation[1] = initialY; initialLocation[2] = initialZ;
        gene = Gene;
        generatecubemasses(initialX, initialY, initialZ);
        generatecubesprings();
    }

    void generatecubemasses(double initialX, double initialY, double initialZ)
    {
        cubemasses.push_back({ mass,{initialX + length / 2,initialY + length / 2,initialZ},{0,0,0},{0,0,0} });
        cubemasses.push_back({ mass,{initialX + length / 2,initialY - length / 2,initialZ},{0,0,0},{0,0,0} });
        cubemasses.push_back({ mass,{initialX - length / 2,initialY - length / 2,initialZ},{0,0,0},{0,0,0} });
        cubemasses.push_back({ mass,{initialX - length / 2,initialY + length / 2,initialZ},{0,0,0},{0,0,0} });
        cubemasses.push_back({ mass,{initialX + length / 2,initialY + length / 2,initialZ + length},{0,0,0},{0,0,0} });
        cubemasses.push_back({ mass,{initialX + length / 2,initialY - length / 2,initialZ + length},{0,0,0},{0,0,0} });
        cubemasses.push_back({ mass,{initialX - length / 2,initialY - length / 2,initialZ + length},{0,0,0},{0,0,0} });
        cubemasses.push_back({ mass,{initialX - length / 2,initialY + length / 2,initialZ + length},{0,0,0},{0,0,0} });
        cubemasses.push_back({ mass,{initialX - 3 * length / 2,initialY + length / 2,initialZ + length},{0,0,0},{0,0,0} });
        cubemasses.push_back({ mass,{initialX - 3 * length / 2,initialY - length / 2,initialZ + length},{0,0,0},{0,0,0} });
        cubemasses.push_back({ mass,{initialX - 3 * length / 2,initialY + length / 2,initialZ},{0,0,0},{0,0,0} });
        cubemasses.push_back({ mass,{initialX - 3 * length / 2,initialY - length / 2,initialZ},{0,0,0},{0,0,0} });
    }

    void generatecubesprings()
    {
        double face_diag = sqrt(pow(length, 2) + pow(length, 2));
        double across_diag = sqrt(pow(length, 2) + pow(length, 2) + pow(length, 2));
      
        cubesprings.push_back({ springConstant,length,   0,1 });
        cubesprings.push_back({ springConstant,face_diag,0,2 });
        cubesprings.push_back({ springConstant,length,   0,3 });
        cubesprings.push_back({ springConstant,length,   0,4 });
        cubesprings.push_back({ springConstant,face_diag,0,5 });
        cubesprings.push_back({ springConstant,across_diag,0,6 });
        cubesprings.push_back({ springConstant,face_diag,0,7 });
        cubesprings.push_back({ springConstant,length,    1,2 });
        cubesprings.push_back({ springConstant,face_diag, 1,3 });
        cubesprings.push_back({ springConstant,face_diag, 1,4 });
        cubesprings.push_back({ springConstant,length,   1,5 });
        cubesprings.push_back({ springConstant,face_diag,1,6 });
        cubesprings.push_back({ springConstant,across_diag,1,7 });
        cubesprings.push_back({ springConstant,length,2,3 });
        cubesprings.push_back({ springConstant,across_diag,2,4 });
        cubesprings.push_back({ springConstant,face_diag,2,5 });
        cubesprings.push_back({ springConstant,length,   2,6 });
        cubesprings.push_back({ springConstant,face_diag,2,7 });
        cubesprings.push_back({ springConstant,face_diag,3,4 });
        cubesprings.push_back({ springConstant,across_diag,3,5 });
        cubesprings.push_back({ springConstant,face_diag,3,6 });
        cubesprings.push_back({ springConstant,length,   3,7 });
        cubesprings.push_back({ springConstant,length,   4,5 });
        cubesprings.push_back({ springConstant,face_diag,4,6 });
        cubesprings.push_back({ springConstant,length,   4,7 });
        cubesprings.push_back({ springConstant,length,   5,6 });
        cubesprings.push_back({ springConstant,face_diag,5,7 });
        cubesprings.push_back({ springConstant,length,6,7 });
        cubesprings.push_back({ springConstant,length,7,8 });
        cubesprings.push_back({ springConstant,length,8,9 });
        cubesprings.push_back({springConstant,length,6,9 });
        cubesprings.push_back({ springConstant,length,8,10 });
        cubesprings.push_back({ springConstant,length,10,11 });
        cubesprings.push_back({ springConstant,length,9,11 });
        cubesprings.push_back({ springConstant,length,3,10 });
        cubesprings.push_back({ springConstant,length,2,11 });
        cubesprings.push_back({ springConstant,face_diag,3,8 });
        cubesprings.push_back({ springConstant,face_diag,7,10 });
        cubesprings.push_back({ springConstant,face_diag,8,11 });
        cubesprings.push_back({ springConstant,face_diag,9,10 });
        cubesprings.push_back({ springConstant,face_diag,2,9 });
        cubesprings.push_back({ springConstant,face_diag,6,11 });
        cubesprings.push_back({ springConstant,face_diag,3,11 });
        cubesprings.push_back({ springConstant,face_diag,2,10 });
        cubesprings.push_back({ springConstant,face_diag,7,9 });
        cubesprings.push_back({ springConstant,face_diag,6,8 });
        cubesprings.push_back({ springConstant,across_diag,2,8 });
        cubesprings.push_back({ springConstant,across_diag,6,10 });
        cubesprings.push_back({ springConstant,across_diag,7,11 });
        cubesprings.push_back({ springConstant,across_diag,3,9 });

    }

    void cubotDraw()
    {
        glColor3f(1, 0, 0);

        GLUquadric* quad;
        quad = gluNewQuadric();
        for (int i = 0; i < (int)cubemasses.size(); i++) {
            glPushMatrix();
            glMultMatrixf(worldRotation);
            glTranslated(cubemasses[i].pos[0], cubemasses[i].pos[1], cubemasses[i].pos[2]);
            gluSphere(quad, length / 20, 10, 10);
            glPopMatrix();
        }

        for (int i = 0; i < (int)cubesprings.size(); i++) {
            glPushMatrix();
            glMultMatrixf(worldRotation);
            glBegin(GL_LINES);
            glVertex3f(GLfloat(cubemasses[cubesprings[i].m1].pos[0]), GLfloat(cubemasses[cubesprings[i].m1].pos[1]), GLfloat(cubemasses[cubesprings[i].m1].pos[2]));
            glVertex3f(GLfloat(cubemasses[cubesprings[i].m2].pos[0]), GLfloat(cubemasses[cubesprings[i].m2].pos[1]), GLfloat(cubemasses[cubesprings[i].m2].pos[2]));
            glEnd();
            glPopMatrix();
        }


       // glColor3f(0.1, 0.3, 0.8);

       // glBegin(GL_QUADS);
       //// glNormal3f(0, 0, 1);
       // glTexCoord2f(0.0f, 0.0f);    glVertex3f(cubemasses[2].pos[0], cubemasses[2].pos[1], cubemasses[2].pos[2]);
       // glTexCoord2f(1.0f, 0.0f);    glVertex3f(cubemasses[3].pos[0], cubemasses[3].pos[1], cubemasses[3].pos[2]);
       // glTexCoord2f(1.0f, 1.0f);    glVertex3f(cubemasses[7].pos[0], cubemasses[7].pos[1], cubemasses[7].pos[2]);
       // glTexCoord2f(0.0f, 1.0f);    glVertex3f(cubemasses[6].pos[0], cubemasses[6].pos[1], cubemasses[6].pos[2]);
       // glEnd();


       // glColor3f(0.2, 0.9, 0.1);

       // glBegin(GL_QUADS);
       //// glNormal3f(0, 0, -1);
       // glTexCoord2f(0.0f, 0.0f);    glVertex3f(cubemasses[0].pos[0], cubemasses[0].pos[1], cubemasses[0].pos[2]);
       // glTexCoord2f(1.0f, 0.0f);    glVertex3f(cubemasses[1].pos[0], cubemasses[1].pos[1], cubemasses[1].pos[2]);
       // glTexCoord2f(1.0f, 1.0f);    glVertex3f(cubemasses[5].pos[0], cubemasses[5].pos[1], cubemasses[5].pos[2]);
       // glTexCoord2f(0.0f, 1.0f);    glVertex3f(cubemasses[4].pos[0], cubemasses[4].pos[1], cubemasses[4].pos[2]);
       // glEnd();


       // glColor3f(0.6, 0.2, 0.3);

       // glBegin(GL_QUADS);
       //// glNormal3f(1, 0, 0);
       // glTexCoord2f(0.0f, 0.0f);    glVertex3f(cubemasses[0].pos[0], cubemasses[0].pos[1], cubemasses[0].pos[2]);
       // glTexCoord2f(1.0f, 0.0f);    glVertex3f(cubemasses[3].pos[0], cubemasses[3].pos[1], cubemasses[3].pos[2]);
       // glTexCoord2f(1.0f, 1.0f);    glVertex3f(cubemasses[7].pos[0], cubemasses[7].pos[1], cubemasses[7].pos[2]);
       // glTexCoord2f(0.0f, 1.0f);    glVertex3f(cubemasses[4].pos[0], cubemasses[4].pos[1], cubemasses[4].pos[2]);
       // glEnd();



       // glColor3f(0.7, 0.0, 0.5);

       // glBegin(GL_QUADS);
       // glNormal3f(-1, 0, 0);
       // glTexCoord2f(0.0f, 0.0f);    glVertex3f(cubemasses[2].pos[0], cubemasses[2].pos[1], cubemasses[2].pos[2]);
       // glTexCoord2f(1.0f, 0.0f);    glVertex3f(cubemasses[1].pos[0], cubemasses[1].pos[1], cubemasses[1].pos[2]);
       // glTexCoord2f(1.0f, 1.0f);    glVertex3f(cubemasses[5].pos[0], cubemasses[5].pos[1], cubemasses[5].pos[2]);
       // glTexCoord2f(0.0f, 1.0f);    glVertex3f(cubemasses[6].pos[0], cubemasses[6].pos[1], cubemasses[6].pos[2]);
       // glEnd();



       // glColor3f(0, 0.4, 0.4);
       // glBegin(GL_QUADS);
       // glNormal3f(0, 1, 0);
       // glTexCoord2f(0.0f, 0.0f);    glVertex3f(cubemasses[4].pos[0], cubemasses[4].pos[1], cubemasses[4].pos[2]);
       // glTexCoord2f(1.0f, 0.0f);    glVertex3f(cubemasses[5].pos[0], cubemasses[5].pos[1], cubemasses[5].pos[2]);
       // glTexCoord2f(1.0f, 1.0f);    glVertex3f(cubemasses[6].pos[0], cubemasses[6].pos[1], cubemasses[6].pos[2]);
       // glTexCoord2f(0.0f, 1.0f);    glVertex3f(cubemasses[7].pos[0], cubemasses[7].pos[1], cubemasses[7].pos[2]);
       // glEnd();


       // glColor3f(0.5, 0.3, 0);

       // glBegin(GL_QUADS);
       // glNormal3f(0, -1, 0);
       // glTexCoord2f(0.0f, 0.0f);    glVertex3f(cubemasses[0].pos[0], cubemasses[0].pos[1], cubemasses[0].pos[2]);
       // glTexCoord2f(1.0f, 0.0f);    glVertex3f(cubemasses[1].pos[0], cubemasses[1].pos[1], cubemasses[1].pos[2]);
       // glTexCoord2f(1.0f, 1.0f);    glVertex3f(cubemasses[2].pos[0], cubemasses[2].pos[1], cubemasses[2].pos[2]);
       // glTexCoord2f(0.0f, 1.0f);    glVertex3f(cubemasses[3].pos[0], cubemasses[3].pos[1], cubemasses[3].pos[2]);
       // glEnd();


       // glColor3f(0.1, 0.3, 0.0);

       // glBegin(GL_QUADS);
       // glNormal3f(0, 0, 1);
       // glTexCoord2f(0.0f, 0.0f);    glVertex3f(cubemasses[8].pos[0], cubemasses[8].pos[1], cubemasses[8].pos[2]);
       // glTexCoord2f(1.0f, 0.0f);    glVertex3f(cubemasses[7].pos[0], cubemasses[7].pos[1], cubemasses[7].pos[2]);
       // glTexCoord2f(1.0f, 1.0f);    glVertex3f(cubemasses[3].pos[0], cubemasses[3].pos[1], cubemasses[3].pos[2]);
       // glTexCoord2f(0.0f, 1.0f);    glVertex3f(cubemasses[10].pos[0], cubemasses[10].pos[1], cubemasses[10].pos[2]);
       // glEnd();


       // glColor3f(0.0, 0.9, 0.1);

       // glBegin(GL_QUADS);
       // glNormal3f(0, 0, 1);
       // glTexCoord2f(0.0f, 0.0f);    glVertex3f(cubemasses[6].pos[0], cubemasses[6].pos[1], cubemasses[6].pos[2]);
       // glTexCoord2f(1.0f, 0.0f);    glVertex3f(cubemasses[8].pos[0], cubemasses[8].pos[1], cubemasses[8].pos[2]);
       // glTexCoord2f(1.0f, 1.0f);    glVertex3f(cubemasses[7].pos[0], cubemasses[7].pos[1], cubemasses[7].pos[2]);
       // glTexCoord2f(0.0f, 1.0f);    glVertex3f(cubemasses[9].pos[0], cubemasses[9].pos[1], cubemasses[9].pos[2]);
       // glEnd();


       // glColor3f(0.4, 0.8, 0.3);

       // glBegin(GL_QUADS);
       // glNormal3f(1, 0, 0);
       // glTexCoord2f(0.0f, 0.0f);    glVertex3f(cubemasses[6].pos[0], cubemasses[6].pos[1], cubemasses[6].pos[2]);
       // glTexCoord2f(1.0f, 0.0f);    glVertex3f(cubemasses[2].pos[0], cubemasses[2].pos[1], cubemasses[2].pos[2]);
       // glTexCoord2f(1.0f, 1.0f);    glVertex3f(cubemasses[9].pos[0], cubemasses[9].pos[1], cubemasses[9].pos[2]);
       // glTexCoord2f(0.0f, 1.0f);    glVertex3f(cubemasses[11].pos[0], cubemasses[11].pos[1], cubemasses[11].pos[2]);
       // glEnd();



       // glColor3f(0.6, 0.1, 0.5);

       // glBegin(GL_QUADS);
       // glNormal3f(-1, 0, 0);
       // glTexCoord2f(0.0f, 0.0f);    glVertex3f(cubemasses[8].pos[0], cubemasses[8].pos[1], cubemasses[8].pos[2]);
       // glTexCoord2f(1.0f, 0.0f);    glVertex3f(cubemasses[9].pos[0], cubemasses[9].pos[1], cubemasses[9].pos[2]);
       // glTexCoord2f(1.0f, 1.0f);    glVertex3f(cubemasses[10].pos[0], cubemasses[10].pos[1], cubemasses[10].pos[2]);
       // glTexCoord2f(0.0f, 1.0f);    glVertex3f(cubemasses[11].pos[0], cubemasses[11].pos[1], cubemasses[11].pos[2]);
       // glEnd();



       // glColor3f(0, 0.3, 0.4);
       // glBegin(GL_QUADS);
       // glNormal3f(0, 1, 0);
       // glTexCoord2f(0.0f, 0.0f);    glVertex3f(cubemasses[10].pos[0], cubemasses[10].pos[1], cubemasses[10].pos[2]);
       // glTexCoord2f(1.0f, 0.0f);    glVertex3f(cubemasses[11].pos[0], cubemasses[11].pos[1], cubemasses[11].pos[2]);
       // glTexCoord2f(1.0f, 1.0f);    glVertex3f(cubemasses[3].pos[0], cubemasses[3].pos[1], cubemasses[3].pos[2]);
       // glTexCoord2f(0.0f, 1.0f);    glVertex3f(cubemasses[2].pos[0], cubemasses[2].pos[1], cubemasses[2].pos[2]);
       // glEnd();




        
        double x = 0; double y = 0; double z = 0;
        for (int j = 1; j < 5; j++) {
            x = x + cubemasses[j].pos[0];
            y = y + cubemasses[j].pos[1];
            z = z + cubemasses[j].pos[2];
        }
        x = x / 4;
        y = y / 4;
        z = z / 4;
        glPushMatrix();
        glMultMatrixf(worldRotation);
        glBegin(GL_LINES);
        glColor3f(0.0, 1.0, 0.0);
        glVertex3f(GLfloat(x), GLfloat(y), GLfloat(z));
        glVertex3f(GLfloat(initialLocation[0]), GLfloat(initialLocation[1]), GLfloat(initialLocation[2]));
        glEnd();
        glPopMatrix();
        Frames++;
        GLint t = glutGet(GLUT_ELAPSED_TIME);
        if (t - T0 >= 1000) {
            GLfloat seconds = (t - T0) / 1000.0;
            fps = Frames / seconds;
            
            T0 = t;
            Frames = 0;
        }
    }

    void cubotUpdate()
    {
      
        std::vector<std::vector<double>> cubotForces((int)cubemasses.size(), std::vector<double>(3));
        if (T > 1.0) {
            cubesprings[cubesprings.size() - 4].L_0 =  gene[0].A  + gene[0].B * cos( PI/ 4 * T + gene[0].C);
            cubesprings[cubesprings.size() - 3].L_0 = gene[1].A  + gene[1].B * cos(PI/ 4* T + gene[1].C);
            cubesprings[cubesprings.size() - 2].L_0 = gene[2].A  + gene[2].B * cos(PI/ 4 * T + gene[2].C);
            cubesprings[cubesprings.size() - 1].L_0 = gene[3].A  + gene[3].B * cos(PI/ 4 * T + gene[3].C);
        }

        
        #pragma omp parallel 
        #pragma omp for
        for (int i = 0; i < (int)cubesprings.size(); i++) {
            Mass mass1 = cubemasses[cubesprings[i].m1];
            Mass mass2 = cubemasses[cubesprings[i].m2];
            double positionDiff[3] = { mass2.pos[0] - mass1.pos[0], mass2.pos[1] - mass1.pos[1], mass2.pos[2] - mass1.pos[2] };
            double L = calcMag(positionDiff, 3);
            double force = cubesprings[i].k * fabs(cubesprings[i].L_0 - L);
            double direction[3] = { positionDiff[0] / L, positionDiff[1] / L, positionDiff[2] / L };
          
            if (L > cubesprings[i].L_0) {
                cubotForces[cubesprings[i].m1][0] = cubotForces[cubesprings[i].m1][0] + direction[0] * force;
                cubotForces[cubesprings[i].m1][1] = cubotForces[cubesprings[i].m1][1] + direction[1] * force;
                cubotForces[cubesprings[i].m1][2] = cubotForces[cubesprings[i].m1][2] + direction[2] * force;
                cubotForces[cubesprings[i].m2][0] = cubotForces[cubesprings[i].m2][0] - direction[0] * force;
                cubotForces[cubesprings[i].m2][1] = cubotForces[cubesprings[i].m2][1] - direction[1] * force;
                cubotForces[cubesprings[i].m2][2] = cubotForces[cubesprings[i].m2][2] - direction[2] * force;
            }
      
            else if (L < cubesprings[i].L_0) {
                cubotForces[cubesprings[i].m1][0] = cubotForces[cubesprings[i].m1][0] - direction[0] * force;
                cubotForces[cubesprings[i].m1][1] = cubotForces[cubesprings[i].m1][1] - direction[1] * force;
                cubotForces[cubesprings[i].m1][2] = cubotForces[cubesprings[i].m1][2] - direction[2] * force;
                cubotForces[cubesprings[i].m2][0] = cubotForces[cubesprings[i].m2][0] + direction[0] * force;
                cubotForces[cubesprings[i].m2][1] = cubotForces[cubesprings[i].m2][1] + direction[1] * force;
                cubotForces[cubesprings[i].m2][2] = cubotForces[cubesprings[i].m2][2] + direction[2] * force;
            }
        }

        #pragma omp parallel 
        #pragma omp for
        for (int i = 0; i < (int)cubemasses.size(); i++) {
           
            cubotForces[i][2] = cubotForces[i][2] - cubemasses[i].m * gravity;
           
            if (cubemasses[i].pos[2] <= 0) {
                cubotForces[i][2] = cubotForces[i][2] + restoreConstant * fabs(cubemasses[i].pos[2]);
               
                double F_h = sqrt(pow(cubotForces[i][0], 2) + pow(cubotForces[i][1], 2));
                double F_v = cubotForces[i][2];
                if (F_h < F_v * frictionCoefficient) {
                    cubotForces[i][0] = 0;
                    cubotForces[i][1] = 0;
                    cubemasses[i].vel[0] = 0;
                    cubemasses[i].vel[1] = 0;
                }
            }
          
            cubemasses[i].acc[0] = cubotForces[i][0] / cubemasses[i].m;
            cubemasses[i].acc[1] = cubotForces[i][1] / cubemasses[i].m;
            cubemasses[i].acc[2] = cubotForces[i][2] / cubemasses[i].m;
        
            cubemasses[i].vel[0] = dampingConstant * (cubemasses[i].vel[0] + cubemasses[i].acc[0] * timeStep);
            cubemasses[i].vel[1] = dampingConstant * (cubemasses[i].vel[1] + cubemasses[i].acc[1] * timeStep);
            cubemasses[i].vel[2] = dampingConstant * (cubemasses[i].vel[2] + cubemasses[i].acc[2] * timeStep);
            
            cubemasses[i].pos[0] = cubemasses[i].pos[0] + cubemasses[i].vel[0] * timeStep;
            cubemasses[i].pos[1] = cubemasses[i].pos[1] + cubemasses[i].vel[1] * timeStep;
            cubemasses[i].pos[2] = cubemasses[i].pos[2] + cubemasses[i].vel[2] * timeStep;
        }
      
        T = T + timeStep;
    }


};

class Simulation {
private:
    int populationSize;
    std::vector<double> populationDistance;
    std::vector<std::vector<Gene>> populationGene;
    std::vector<std::vector<Gene>> newPopulationGene;
    std::vector<cubot> cubots;
public:
    double averageDistance;
    double maxDistance;

    Simulation(int popSize)
    {
        populationSize = popSize;
        generateGenes();
     
        generatecubots();
        popDis.open("populationDistance.txt"); 
        bestGene.open("bestGene.txt"); 
        
    }

    void startSim(double time) {
        if (T < time) {
            simUpdate();
            simDraw();
        }
        else {
            calculatePopulationDistance();
            selection();
            crossOver();
            populationGene.clear(); populationGene.shrink_to_fit();
            populationGene = newPopulationGene;
            generationNumber++;
            cubots.clear(); cubots.shrink_to_fit();
            generatecubots();
            T = 0;
        }
    }

    void selection() {
        std::vector<int> index = sort(populationDistance);
        newPopulationGene.clear();
        newPopulationGene.shrink_to_fit();
        for (int i = 0; i < index.size() / 2; i++) {
           
            newPopulationGene.push_back(populationGene[index[i]]);
        }
        for (int i = 0; i < newPopulationGene[0].size(); i++) {
            bestGene << newPopulationGene[0][i].A << " " << newPopulationGene[0][i].B << " " << newPopulationGene[0][i].C << " ";
        }
        bestGene << "\n";
    }

    void crossOver() {
        #pragma omp parallel 
        #pragma omp for
        for (int n = 0; n < populationGene.size() / 2; n++) {
            int parentIndex1 = rand() % static_cast<int>(newPopulationGene.size());
            int parentIndex2 = rand() % static_cast<int>(newPopulationGene.size());
            std::vector<double> parent1, parent2;
            for (int i = 0; i < newPopulationGene[parentIndex1].size(); i++) {
                parent1.push_back(newPopulationGene[parentIndex1][i].A);
                parent1.push_back(newPopulationGene[parentIndex1][i].B);
                parent1.push_back(newPopulationGene[parentIndex1][i].C);
                parent2.push_back(newPopulationGene[parentIndex2][i].A);
                parent2.push_back(newPopulationGene[parentIndex2][i].B);
                parent2.push_back(newPopulationGene[parentIndex2][i].C);
            }
            int crossOverPoint1 = rand() % static_cast<int>(parent1.size());
            int crossOverPoint2 = rand() % static_cast<int>(parent1.size());
            if (crossOverPoint2 < crossOverPoint1) {
                int temp = crossOverPoint1;
                crossOverPoint1 = crossOverPoint2;
                crossOverPoint2 = temp;
            }
            std::vector<double> offSpring1, offSpring2;
            for (int i = 0; i < crossOverPoint1; i++) {
                offSpring1.push_back(parent1[i]);
                offSpring2.push_back(parent2[i]);
            }
            for (int i = crossOverPoint1; i < crossOverPoint2; i++) {
                offSpring1.push_back(parent2[i]);
                offSpring2.push_back(parent1[i]);
            }
            for (int i = crossOverPoint2; i < parent1.size(); i++) {
                offSpring1.push_back(parent1[i]);
                offSpring2.push_back(parent2[i]);
            }
            offSpring1 = mutation(offSpring1);
            offSpring2 = mutation(offSpring2);
            std::vector<Gene> offSpringGene1, offSpringGene2;
            Gene temp1, temp2;
            for (int i = 0; i < offSpring1.size(); i = i + 3) {
                temp1.A = offSpring1[i];
                temp1.B = offSpring1[i + 1];
                temp1.C = offSpring1[i + 2];
                temp2.A = offSpring2[i];
                temp2.B = offSpring2[i + 1];
                temp2.C = offSpring2[i + 2];
                offSpringGene1.push_back(temp1);
                offSpringGene2.push_back(temp2);
            }
            newPopulationGene.push_back(offSpringGene1);
            newPopulationGene.push_back(offSpringGene2);
        }
    }

    std::vector<double> mutation(std::vector<double> offSpring) {
        for (int i = 0; i < offSpring.size(); i++) {
            double mutationProbability = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 1.0));
            if (mutationProbability > 0.9) {
                offSpring[i] = -2 * PI + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (4 * PI)));
            }
        }
        return offSpring;

    }

    void generateBestGene() {
        for (int i = 0; i < populationSize; i++) {
            std::vector<Gene> tempVec1;
            Gene temp1{ 7.25, -6.729, 2.43 };
            tempVec1.push_back(temp1);
            Gene temp2{ -0.0439628, -0.281218, 1.16934 };
            tempVec1.push_back(temp2);
            Gene temp3{ -0.00115234, -2.94204, 6.13413 };
            tempVec1.push_back(temp3);
            Gene temp4{ -5.13257, -0.0672038, 1.83822 };
            tempVec1.push_back(temp4);
            populationGene.push_back(tempVec1);
        }
    }

    void generateGenes() {
        srand(time(0));
        for (int i = 0; i < populationSize; i++) {
            std::vector<Gene> tempVec;
            for (int j = 0; j < 4; j++) {
                double A = -2 * PI + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (4 * PI))); 
                double B = -2 * PI + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (4 * PI))); 
                double C = -2 * PI + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (4 * PI))); 
                Gene temp{ A,B,C };
                tempVec.push_back(temp);
            }
            populationGene.push_back(tempVec);
        }
    }

    void generatecubots() {
        #pragma omp parallel 
        #pragma omp for
        for (int i = 0; i < populationSize; i++) {
            double X = -20 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 40));
            double Y = -20 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 40));
            
            cubots.push_back(cubot(X, Y, 7.0, populationGene[i]));
        }
    }

    void simUpdate() {
        #pragma omp parallel 
        #pragma omp for

        for (int i = 0; i < populationSize; i++) {
            cubots[i].cubotUpdate();
        }
    }

    void simDraw() {
        for (int i = 0; i < populationSize; i++) {
            cubots[i].cubotDraw();
        }
    }

    void calculatePopulationDistance() {
        populationDistance.clear();
        populationDistance.shrink_to_fit();
        for (int i = 0; i < populationSize; i++) {
            double x = 0; double y = 0;
            for (int j = 1; j < 5; j++) {
                x = x + cubots[i].cubemasses[j].pos[0];
                y = y + cubots[i].cubemasses[j].pos[1];
            }
            x = x / 4;
            y = y / 4;
            double distance[2] = { fabs(x - cubots[i].initialLocation[0]), fabs(y - cubots[i].initialLocation[1]) };
            double distancecalcMag = calcMag(distance, 2);
            populationDistance.push_back(distancecalcMag);
        }
        averageDistance = 0;
        maxDistance = 0;
        std::cout << "##################" << std::endl;
        for (int i = 0; i < populationSize; i++) {
            averageDistance = averageDistance + populationDistance[i];
            maxDistance = std::max(maxDistance, populationDistance[i]);
           
            popDis << populationDistance[i] << " ";
        }
        popDis << "\n";
        averageDistance = averageDistance / populationSize;
        std::cout << "Maximum Distance: " << maxDistance << std::endl;
        std::cout << "Average Distance: " << averageDistance << std::endl;
    }
};

Simulation sim1(cubotNumber);

void Print(const char* format, ...)
{
    char    buf[LEN];
    char* ch = buf;
    va_list args;
   
    va_start(args, format);
    vsnprintf(buf, LEN, format, args);
    va_end(args);
    
    while (*ch)
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, *ch++);
}



void drawGrid() {

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);



    glRotatef(0, 1.0, 0.0, 2.0);

    glRotatef(0, 0.0, 1.0, 2.0);

    glColor3f(0.0, 1.0, 2.0);
    glBegin(GL_LINES);
    for (GLfloat i = -30; i < 30; i += 1)
    {
        glVertex3f(i, 0, 30);
        glVertex3f(i, 0, -30);
        glVertex3f(30, 0, i);
        glVertex3f(-30, 0, i);
    }
    glEnd();


}


void display()
{
    const double len = 2.0; 
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 
    glEnable(GL_DEPTH_TEST);
    
    glLoadIdentity();
   
    double Ex = -2 * dim * Sin(th) * Cos(ph);
    double Ey = +2 * dim * Sin(ph);
    double Ez = +2 * dim * Cos(th) * Cos(ph);
    gluLookAt(Ex, Ey, Ez, 0, 0, 0, 0, Cos(ph), 0);

   

    drawGrid();
    glColor3f(1, 1, 1);

 
    sim1.startSim(simulationTime);

    glRasterPos3d(0.0, 2, 0.0);
   

    glColor3f(1, 1, 1);

    glFlush();

    glutSwapBuffers();
}


void special(int key, int x, int y)
{
    
    if (key == GLUT_KEY_RIGHT)
        th += 5;

    else if (key == GLUT_KEY_LEFT)
        th -= 5;
    else if (key == GLUT_KEY_UP)
    {
        if (ph + 5 < 90)
        {
            ph += 5;
        }
    }
  
    else if (key == GLUT_KEY_DOWN)
    {
        if (ph - 5 > 0)
        {
            ph -= 5;
        }
    }
 
    th %= 360;
    ph %= 360;

    glutPostRedisplay();
}


void Project(double fov, double asp, double dim)
{

    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();

    if (fov)
        gluPerspective(fov, asp, dim / 16, 16 * dim);

    else
        glOrtho(-asp * dim, asp * dim, -dim, +dim, -dim, +dim);

    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();
}

void reshape(int width, int height)
{

    asp = (height > 0) ? (double)width / height : 1;

    glViewport(0, 0, width, height);

    Project(fov, asp, dim);
}

void idle()
{
    glutPostRedisplay();
}

int main(int argc, char* argv[]) {
   
    glutInit(&argc, argv);
 
    glutInitWindowSize(1000, 800);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
    glutCreateWindow("cubot_sim_jb4512");
    glutIdleFunc(idle);
    glutDisplayFunc(display);

    glutReshapeFunc(reshape);
    glutSpecialFunc(special);

    glutMainLoop();

    return 0;
}