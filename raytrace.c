#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

//#define DEBUG
#define AUTHOR "CBLAZER"
#define RGB_NUMBER 255
#define RECURSION_LIMIT 7
#define PI 3.14159


//struct for holding objects
typedef struct Object{
    char *type;
    double *difColor;
    double *specColor;
    double *position;
    double *normal;
    double radius;
    double reflect;
    double refract;
    double ior;
} Object;

Object objects[128]; //128 objects per assignment

//struct for holding lights
typedef struct Light{
    char *type;
    double *position;
    double *color;
    double radialA0;
    double radialA1;
    double radialA2;
    double theta;
    double angularA0;
    double *direction;
} Light;

Light lights[128]; //array of lights


typedef struct Pixel{
    double red;
    double green;
    double blue;
} Pixel;



int line = 1;
int Width;
int Height;
int color;
double cameraWidth;
double cameraHeight;
Pixel *viewPlane; //2d array to store pixel color data from raycast
double black[3] = {0,0,0};

//number maintenance
int next_c(FILE* json) {
  int c = fgetc(json);
#ifdef DEBUG
  printf("next_c: '%c'\n", c);
#endif
  if (c == '\n') {
    line += 1;
  }
  if (c == EOF) {
    fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", line);
    exit(1);
  }
  return c;
}

//expect_c() checks that the next character is d.  If it is not it emits an error.
void expect_c(FILE* json, int d) {
  int c = next_c(json);
  if (c == d) return;
  fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, line);
  exit(1);
}

//skip_ws() skips white space in the file.
void skip_ws(FILE* json) {
  int c = next_c(json);
  while (isspace(c)) {
    c = next_c(json);
  }
  ungetc(c, json);
}

//next_string() gets the next string from the file handle and emits an error if a string can not be obtained.
char* next_string(FILE* json) {
  char buffer[129];
  int c = next_c(json);
  if (c != '"') {
    fprintf(stderr, "Error: Expected string on line %d.\n", line);
    exit(1);
  }
  c = next_c(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
      exit(1);
    }
    if (c == '\\') {
      fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
      exit(1);
    }
    if (c < 32 || c > 126) {
      fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
      exit(1);
    }
    buffer[i] = c;
    i += 1;
    c = next_c(json);
  }
  buffer[i] = 0;
  return strdup(buffer);
}

double next_number(FILE* json) {
  double value;
  fscanf(json, "%lf", &value);
  return value;
}

double* next_vector(FILE* json) {
  double* v = malloc(3*sizeof(double));
  expect_c(json, '[');
  skip_ws(json);
  v[0] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[1] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[2] = next_number(json);
  skip_ws(json);
  expect_c(json, ']');
  return v;
}

void read_scene(char* filename) {
  int c;
  FILE* json = fopen(filename, "r");

  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", filename);
    exit(1);
  }

  skip_ws(json);

  //Find the beginning of the list
  expect_c(json, '[');

  skip_ws(json);

  //Find the objects
  int objectIndex = -1;
  int lightIndex = -1;
  int lightPosition = 0; //variable to keep track of light or object because they both have a "position" variable
  while (1) {
    c = fgetc(json);
    if (c == ']') {
      fprintf(stderr, "Error: This is the worst scene file EVER.\n");
      fclose(json);
      return;
    }
    if (c == '{') {
      skip_ws(json);

      //Parse the object
      char* key = next_string(json);
      if (strcmp(key, "type") != 0) {
        fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", line);
        exit(1);
      }

      skip_ws(json);

      expect_c(json, ':');

      skip_ws(json);

      char* value = next_string(json);

      if (strcmp(value, "camera") == 0) {
        //Do nothing, camera isn't an object in the scene.
      } else if (strcmp(value, "sphere") == 0) {
        objectIndex++;
        objects[objectIndex].type = "sphere";
        lightPosition = 0;
      } else if (strcmp(value, "plane") == 0) {
        objectIndex++;
        objects[objectIndex].type = "plane";
        lightPosition = 0;
      } else if (strcmp(value, "light") == 0) {
        lightIndex++;
        lights[lightIndex].type = "light";
        lightPosition = 1;
      } else {
        fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
        exit(1);
      }

      skip_ws(json);

      while (1) {

        c = next_c(json);
        if (c == '}') {
          //stop parsing this object
          break;
        } else if (c == ',') {
          //read another field
          skip_ws(json);
          char* key = next_string(json);
          skip_ws(json);
          expect_c(json, ':');
          skip_ws(json);
          if (strcmp(key, "width") == 0) {
            cameraWidth = next_number(json);
          }
          else if (strcmp(key, "height") == 0) {
            cameraHeight = next_number(json);
          }
          else if (strcmp(key, "radius") == 0) {
            objects[objectIndex].radius = next_number(json);
          }
          else if (strcmp(key, "reflectivity") == 0) {
            objects[objectIndex].reflect = next_number(json);
          }
          else if (strcmp(key, "refractivity") == 0) {
            objects[objectIndex].refract = next_number(json);
          }
          else if (strcmp(key, "ior") == 0) {
            objects[objectIndex].ior = next_number(json);
          }
          else if (strcmp(key, "diffuse_color") == 0) {
            objects[objectIndex].difColor = next_vector(json);
          }
          else if (strcmp(key, "specular_color") == 0) {
            objects[objectIndex].specColor = next_vector(json);
          }
          else if (strcmp(key, "position") == 0) {
            if(lightPosition == 1){
              lights[lightIndex].position = next_vector(json);
            }
            else{
              objects[objectIndex].position = next_vector(json);
            }
          }
          else if (strcmp(key, "normal") == 0) {
            objects[objectIndex].normal = next_vector(json);
          }
          else if (strcmp(key, "color") == 0) {
            lights[lightIndex].color = next_vector(json);
          }
          else if (strcmp(key, "radial-a0") == 0) {
            lights[lightIndex].radialA0 = next_number(json);
          }
          else if (strcmp(key, "radial-a1") == 0) {
            lights[lightIndex].radialA1 = next_number(json);
          }
          else if (strcmp(key, "radial-a2") == 0) {
            lights[lightIndex].radialA2 = next_number(json);
          }
          else if (strcmp(key, "theta") == 0) {
            lights[lightIndex].theta = next_number(json);
          }
          else if (strcmp(key, "angular-a0") == 0) {
            lights[lightIndex].angularA0 = next_number(json);
          }
          else if (strcmp(key, "direction") == 0) {
            lights[lightIndex].direction = next_vector(json);
          }
          else {
            fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n",
              key, line);
          }
          skip_ws(json);
        } else {
          fprintf(stderr, "Error: Unexpected value on line %d\n", line);
          exit(1);
        }
      }
      skip_ws(json);
      c = next_c(json);
      if (c == ',') {

  skip_ws(json);
      } else if (c == ']') {
  fclose(json);
  return;
      } else {
        fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
        exit(1);
      }
    }
  }
}

// added my own sqr function like c's built in sqrt
double sqr(double num){
  return num*num;
}

double tClosestApproachPlane(double *normal, double *origin, double *position, double *lookUVector) {
  double out = -(normal[0] * (origin[0] - position[0]) + normal[1] * (origin[1] - position[1]) + normal[2] * (origin[2] - position[2])) / (normal[0] * lookUVector[0] + normal[1] * lookUVector[1] + normal[2] * lookUVector[2]);
  return(out);
}

double tClosestApproachSphere(double *in1, double *in2){
  double out = (in1[0] * in2[0] + in1[1] * in2[1] + in1[2] * in2[2]) / (sqr(in1[0]) + sqr(in1[1]) + sqr(in1[2]));
  return(out);
}

void vectorAdd(double* in1, double* in2, double* out){
  out[0] = in1[0] + in2[0];
  out[1] = in1[1] + in2[1];
  out[2] = in1[2] + in2[2];
}

void vectorSub(double* in1, double* in2, double* out){
  out[0] = in1[0] - in2[0];
  out[1] = in1[1] - in2[1];
  out[2] = in1[2] - in2[2];
}

void vectorMult(double* in1, double in2, double* out){

  out[0] = in1[0] * in2;
  out[1] = in1[1] * in2;
  out[2] = in1[2] * in2;
}

void vectorDiv(double* in1, double in2, double* out){
  out[0] = in1[0] / in2;
  out[1] = in1[1] / in2;
  out[2] = in1[2] / in2;
}

double vectorDot(double* in1, double* in2){
  double out = (in1[0] * in2[0] + in1[1] * in2[1] + in1[2] * in2[2]);
  return(out);
}

double vectorMag(double *vector) {
  return sqrt(sqr(vector[0]) + sqr(vector[1]) + sqr(vector[2]));
}

void vectorUnit(double *in, double *out){
  vectorDiv(in, vectorMag(in), out);
}


double distance(double* in1, double* in2){
  double out = sqrt(pow(in1[0] - in2[0], 2) + pow(in1[1] - in2[1], 2) + pow(in1[2] - in2[2], 2));
  return(out);
}

void vectorReflect(double* normal, double* lightV, double* reflect){

  vectorUnit(normal, normal);
  vectorUnit(lightV, lightV);
  double temp[3] = {0, 0, 0};
  double dotProduct = vectorDot(normal, lightV);
  vectorMult(lightV, dotProduct, temp);
  vectorMult(lightV, 2.0, lightV);
  vectorSub(normal, lightV, reflect);

}

double planeIntersect(double* Ro, double* Rd, double* position, double* normal){

    vectorUnit(normal, normal);
    double a = normal[0];
    double b = normal[1];
    double c = normal[2];
    vectorUnit(Rd, Rd);

    double x0 = position[0];
    double y0 = position[1];
    double z0 = position[2];

    double d = -(a*x0 + b*y0 + c*z0);

    double e = (a*Rd[0] + b*Rd[1] + c*Rd[2]);
    if(e == 0.0) return -1;
    double t = -(a*Ro[0] + b*Ro[1] + c*Ro[2] + d)/(a*Rd[0] + b*Rd[1] + c*Rd[2]);

    return t;
}

double sphereIntersect(double* Ro, double* Rd, double* position, double radius){

    double x = position[0];
    double y = position[1];
    double z = position[2];
    vectorUnit(Rd, Rd);

    double a = sqr(Rd[0])+sqr(Rd[1])+sqr(Rd[2]);
    double b = 2*(Rd[0]*(Ro[0]-x) + Rd[1]*(Ro[1]-y) + Rd[2]*(Ro[2]-z));
    double c = (sqr((Ro[0]-x)) + sqr((Ro[1]-y)) + sqr((Ro[2]-z)) - sqr(radius));

    double d = (sqr(b) - 4*a*c);

    //no intersection
    if(d < 0.0){
      //return INFINITY;
      return INFINITY;
    }

    double t = ((-b - sqrt(sqr(b) - 4.0*c*a))/(2.0*a));
    if(t > 0){
      return t;
    }
    else{
      t = ((-b + sqrt(sqr(b) - 4.0*c*a))/(2.0*a));
      return t;
    }

    return -1;

}

double fAng(Light a, double* newRd){
    double* lightDirection = malloc(sizeof(double)*3);
    lightDirection[0] = -1*newRd[0];
    lightDirection[1] = -1*newRd[1];
    lightDirection[2] = -1*newRd[2];

    if(a.theta == 0){
      return 1;
    }

    double temp = ((a.direction[0]*lightDirection[0]) +(a.direction[1]*lightDirection[1]) + (a.direction[2]*lightDirection[2]));

    if(temp < (cos(a.theta * PI/180))) return 0;

    free(lightDirection);
    return pow(temp, a.angularA0);
}

double fRad(Light a, double* newRo){
    double* lightV = malloc(sizeof(double)*3);
    vectorSub(a.position, newRo, lightV);
    double b = a.radialA2*sqr(vectorMag(lightV)) + a.radialA1*vectorMag(lightV) + a.radialA0;
    double c = 1*sqr(vectorMag(lightV));
    free(lightV);

    if(b == 0.0){
      return 1/c;
    }
    else{
      return 1/b;
    }
}
double difContribution(int index, double* difColor, Light currentLight, double NdL){
    if(NdL <= 0){
      return 0;
    }
    return difColor[index]*currentLight.color[index]*NdL;
}

double specContribution(int index, double* specColor, Light currentLight, double NdL, double VdR, double ns){
    if(VdR <= 0){
      return 0;
    }
    if(NdL <= 0){
      return 0;
    }
    return specColor[index]*currentLight.color[index]*pow(VdR,ns);
}


double shootD(double *Ro, double *Rd){
    vectorUnit(Rd, Rd);
    double min = INFINITY;

    int objectIndex = 0;
    while(objects[objectIndex].difColor != NULL){
        double t = 0;
        if(strcmp(objects[objectIndex].type, "sphere") == 0){
          t = sphereIntersect(Ro, Rd, objects[objectIndex].position, objects[objectIndex].radius);
        }
        if(strcmp(objects[objectIndex].type, "plane") == 0){
          t = planeIntersect(Ro, Rd, objects[objectIndex].position, objects[objectIndex].normal);
        }
        if(t > 0 && t < min){
          min = t;
        }
        objectIndex++;
    }
    return min;
}
int shootI(double* Ro, double* Rd){
    vectorUnit(Rd, Rd);
    double min = INFINITY;
    int closestObjectIndex = -1;
    int objectIndex = 0;
    while(objects[objectIndex].difColor != NULL){
        double t = 0;
        if(strcmp(objects[objectIndex].type, "sphere") == 0){
                t = sphereIntersect(Ro, Rd, objects[objectIndex].position, objects[objectIndex].radius);
        }
        if(strcmp(objects[objectIndex].type, "plane") == 0){
                t = planeIntersect(Ro, Rd, objects[objectIndex].position, objects[objectIndex].normal);
        }
        if(t > 0 && t < min){
            min = t;
            closestObjectIndex = objectIndex;
        }
        objectIndex++;
    }
    //printf("closestObjectIndex: %d\n", closestObjectIndex);
    return closestObjectIndex;
}

double* shade(int recursionLimit, double *color, int closestObjectIndex, double *Ro, double *Rd, double min, double initIor){

    // stop recursing
    if(recursionLimit == 0 || closestObjectIndex < 0 || min == INFINITY){
        //printf("test\n");
        return black;
    }


    // new value for Ro, instead of camera is is intersection point
    double newRo[3] = {0,0,0};
    newRo[0] = (min * Rd[0]) + Ro[0];
    newRo[1] = (min * Rd[1]) + Ro[1];
    newRo[2] = (min * Rd[2]) + Ro[2];


    int lightIndex = 0;

    while(lights[lightIndex].color != NULL){


        double newRd[3] = {0,0,0};
        vectorSub(lights[lightIndex].position, newRo, newRd);
        vectorUnit(newRd, newRd);

        double distanceToPix = sqrt(sqr(lights[lightIndex].position[0]-newRo[0])+ sqr(lights[lightIndex].position[1]-newRo[1])+ sqr(lights[lightIndex].position[2]-newRo[2]));


        double shadowMin = INFINITY;
        int closestShadowIndex = -1;

        int objectIndex = -1;
        while(objects[objectIndex + 1].difColor != NULL){

          objectIndex++;

          double tShadow = 0;

          if(objectIndex == closestObjectIndex){
            continue;
          }
          if(strcmp(objects[objectIndex].type, "sphere") == 0){
            tShadow = sphereIntersect(newRo, newRd, objects[objectIndex].position, objects[objectIndex].radius);
          }
          if(strcmp(objects[objectIndex].type, "plane") == 0){
            tShadow = planeIntersect(newRo, newRd, objects[objectIndex].position, objects[objectIndex].normal);
          }
          if(min > distanceToPix){
            continue;
          }
          if(tShadow > 0.0 && tShadow < shadowMin){
            shadowMin = tShadow;
            closestShadowIndex = objectIndex;
          }
        }

        if(closestShadowIndex < 0){

          double closeNorm[3] = {0,0,0};

          if(strcmp(objects[closestObjectIndex].type, "sphere") == 0){
            vectorSub(newRo, objects[closestObjectIndex].position, closeNorm);
          }
          if(strcmp(objects[closestObjectIndex].type, "plane") == 0){
            closeNorm[0] = objects[closestObjectIndex].normal[0];
            closeNorm[1] = objects[closestObjectIndex].normal[1];
            closeNorm[2] = objects[closestObjectIndex].normal[2];
          }

          vectorUnit(closeNorm, closeNorm);
          double *vDirection = newRd;
          double reflectionVector[3] = {0,0,0};

          vectorReflect(closeNorm, vDirection, reflectionVector);

          vectorUnit(reflectionVector, reflectionVector);

          double cVector[3] = {0,0,0};
          cVector[0] = Rd[0];
          cVector[1] = Rd[1];
          cVector[2] = Rd[2];
          vectorUnit(cVector, cVector);

          double* currentDif = objects[closestObjectIndex].difColor;
          double* currentSpec = objects[closestObjectIndex].specColor;
          double ns = 20;

          vectorUnit(closeNorm, closeNorm);
          vectorUnit(vDirection, vDirection);


          double VdR = vectorDot(cVector, reflectionVector);
          double NdL = vectorDot(closeNorm, vDirection);

          // reflection
          double rcRo[3] = {0,0,0};
          memcpy(rcRo,newRo, sizeof(double)*3);

          double rcRd[3] = {0,0,0};
          double c1;
          c1 = vectorDot(closeNorm, Rd);
          c1 = -1*c1;
          vectorMult(closeNorm, c1, rcRd);
          vectorMult(rcRd, 2, rcRd);
          vectorAdd(Rd, rcRd, rcRd);
          vectorUnit(rcRd, rcRd);

          rcRo[0] = newRo[0] + rcRd[0]*0.01;
          rcRo[1] = newRo[1] + rcRd[1]*0.01;
          rcRo[2] = newRo[2] + rcRd[2]*0.01;

          double rcMin = shootD(rcRo,rcRd);
          int rcIndex = shootI(rcRo,rcRd);


          // refraction
          double rcRo2[3] = {0,0,0};
          memcpy(rcRo, newRo, sizeof(double)*3);

          double oIor = initIor;
          double newIor = objects[closestObjectIndex].ior;
          double ior = oIor/newIor;


          double c2 = sqrt(1 - sqr(ior) * (1 - sqr(c1)));
          double rcRd2[3] = {0,0,0};
          double temp[3] = {0,0,0};
          vectorMult(Rd, ior, temp);
          vectorMult(closeNorm, (ior*c1-c2), rcRd2);
          vectorAdd(temp,rcRd2,rcRd2);
          vectorMult(rcRd2, -1, rcRd2);
          vectorUnit(rcRd2, rcRd2);

          rcRo2[0] = newRo[0] + rcRd2[0]*0.01;
          rcRo2[1] = newRo[1] + rcRd2[1]*0.01;
          rcRo2[2] = newRo[2] + rcRd2[2]*0.01;

          //printf("Test4\n");

          double rcMin2 = shootD(rcRo2,rcRd2);
          int rcIndex2 = shootI(rcRo2,rcRd2);

          double reflectColor[3] = {0,0,0};
          double refractColor[3] = {0,0,0};

          // printf("recursionLimit: %d\n", recursionLimit);
          // printf("reflectColor: %f, %f, %f\n", reflectColor[0], reflectColor[1], reflectColor[2]);
          //printf("rcIndex1, 2: %d, %d\n", rcIndex, rcIndex2);
          // printf("rcRo: %f, %f, %f\n", rcRo[0], rcRo[1], rcRo[2]);
          // printf("rcRd: %f, %f, %f\n", rcRd[0], rcRd[1], rcRd[2]);
          // printf("rcMin: %f\n", rcMin);
          // printf("ior: %f\n", ior);
          if(objects[closestObjectIndex].reflect != 0){
            shade(recursionLimit-1, reflectColor, rcIndex, rcRo, rcRd, rcMin, ior);
          }
          if(objects[closestObjectIndex].refract != 0){
            shade(recursionLimit-1, refractColor, rcIndex2, rcRo2, rcRd2, rcMin2,ior);
          }
          color[0] += (1 - objects[closestObjectIndex].reflect - objects[closestObjectIndex].refract) * fRad(lights[lightIndex], newRo) * fAng(lights[lightIndex], newRd) * (difContribution(0, currentDif, lights[lightIndex], NdL) + specContribution(0,currentSpec, lights[lightIndex], NdL, VdR, ns)) + ((objects[closestObjectIndex].reflect) * reflectColor[0]) + (objects[closestObjectIndex].refract) * refractColor[0];

          //printf("%f\n", color[0]);
          color[1] += (1 - objects[closestObjectIndex].reflect - objects[closestObjectIndex].refract) * fRad(lights[lightIndex], newRo) * fAng(lights[lightIndex], newRd) * (difContribution(1, currentDif, lights[lightIndex], NdL) + specContribution(1,currentSpec, lights[lightIndex], NdL, VdR, ns)) + ((objects[closestObjectIndex].reflect) * reflectColor[1]) + (objects[closestObjectIndex].refract) * refractColor[1];

          color[2] += (1 - objects[closestObjectIndex].reflect - objects[closestObjectIndex].refract) * fRad(lights[lightIndex], newRo) * fAng(lights[lightIndex], newRd) * (difContribution(2, currentDif, lights[lightIndex], NdL) + specContribution(2,currentSpec, lights[lightIndex], NdL, VdR, ns)) + ((objects[closestObjectIndex].reflect) * reflectColor[2]) + (objects[closestObjectIndex].refract) * refractColor[2];
        }
        lightIndex++;
    }
    return color;
}



void raycast() {

  int i;
  int j;
  int index = 0;
  int objectIndex;

  // camera position
  double Ro[3] = {0.0, 0.0, 0.0};

  //loop through all pixels
  for(i = Width; i > 0; i--){
    for(j = 0; j < Height; j++){

      double x, y, z = 1; //z is always 1 because the view plane is 1 unit away from camera

      //find vector
      x = 0 - (cameraWidth/2) + ((cameraWidth/Width)*(j + 0.5));
      y = 0 - (cameraHeight/2) + ((cameraHeight/Height)*(i + 0.5));

      double Rd[3] = {x, y, z};
      vectorUnit(Rd, Rd);

      double min = INFINITY; //set min so that close objects display over further ones

      double t;
      int closestObjectIndex = -1;

      //loop through all objects
      objectIndex = -1;
      while(objects[objectIndex + 1].difColor != NULL){
        objectIndex++;

        if(strcmp(objects[objectIndex].type, "sphere") == 0){

          t = sphereIntersect(Ro, Rd, objects[objectIndex].position, objects[objectIndex].radius);

            //set new min so that close spheres display over further ones
            if (t > 0 && min >= t){
              min = t;
              closestObjectIndex = objectIndex;
            }
        }
        else if(strcmp(objects[objectIndex].type, "plane") == 0){

          t = planeIntersect(Ro, Rd, objects[objectIndex].position, objects[objectIndex].normal);
          if (t > 0 && min >= t) {
            min = t;
            closestObjectIndex = objectIndex;
          }
        }

      }

      double color[3] = {0,0,0};

      color[0] = 0;
      color[1] = 0;
      color[2] = 0;

      double colorActual[3] = {0,0,0};

      if(min > 0 && min != INFINITY){
        // here is the recursion, gets the color to set the pixel as
        memcpy(colorActual, shade(RECURSION_LIMIT, color, closestObjectIndex, Ro, Rd, min, 1), sizeof(double)*3);
      }

      viewPlane[index].red = colorActual[0];
      viewPlane[index].green = colorActual[1];
      viewPlane[index].blue = colorActual[2];


      index++;
    }
  }
}



void write_scene(char *filename, int format) {
  FILE *ppm = fopen(filename, "wb");
  if (!ppm) {
    fprintf(stderr, "Error opening ppm file\n");
    exit(1);
  }
  //header
  if (format == 6) {
    fprintf(ppm, "P6\n");
  }
  else if (format == 3) {
    fprintf(ppm, "P3\n");
  }
  fprintf(ppm, "# Created by %s\n", AUTHOR);
  fprintf(ppm, "%d %d\n", Width, Height);
  fprintf(ppm, "%d\n", RGB_NUMBER);

  //image data
  int index;
  if (format == 6) {
    for (index = 0; index < Width * Height; index++) {
      color = (int) (viewPlane[index].red * 255); //cast color to int from the double stroed in viewPane
      fwrite(&color, 1, 1, ppm); //red
      color = (int) (viewPlane[index].green * 255);
      fwrite(&color, 1, 1, ppm); //green
      color = (int) (viewPlane[index].blue * 255);
      fwrite(&color, 1, 1, ppm); //blue
    }
  }
  else if (format == 3) {
    for (index = 0; index < Width * Height; index++) {

      color = (int) (viewPlane[index].red * 255);
      fprintf(ppm, "%d\n", color);
      //printf("%d\n", color);
      color = (int) (viewPlane[index].green * 255);
      fprintf(ppm, "%d\n", color);
      //printf("%d\n", color);
      color = (int) (viewPlane[index].blue * 255);
      fprintf(ppm, "%d\n", color);
      //printf("%d\n", color);
    }
  }

  fclose(ppm);
}

void testPrint(){

  printf("testing\n");

  int index = 0;

  while(objects[index].type != NULL){

    if (objects[index].type != NULL) {
      printf("type: %s\n", objects[index].type);
    }

    if (objects[index].radius != 0) {
      printf("radius: %f\n", objects[index].radius);
    }

    if (objects[index].ior != 0) {
      printf("ior: %f\n", objects[index].ior);
    }

    if (objects[index].difColor != NULL) {
      printf("difColor: %f, %f, %f\n", objects[index].difColor[0], objects[index].difColor[1], objects[index].difColor[2]);
    }

    if (objects[index].specColor != NULL) {
      printf("specColor: %f, %f, %f\n", objects[index].specColor[0], objects[index].specColor[1], objects[index].specColor[2]);
    }

    if (objects[index].position != NULL) {
      printf("postion: %f, %f, %f\n", objects[index].position[0], objects[index].position[1], objects[index].position[2]);
    }
    if (objects[index].normal != NULL) {
      printf("normal: %f, %f, %f\n", objects[index].normal[0], objects[index].normal[1], objects[index].normal[2]);
    }
    index++;
    printf("\n\n");
  }


  index = 0;

  while (lights[index].type != NULL) {

    if (lights[index].type != NULL) {
      printf("type: %s\n", lights[index].type);
    }
    if (lights[index].position != NULL) {
      printf("postion: %f, %f, %f\n", lights[index].position[0], lights[index].position[1], lights[index].position[2]);
    }
    if (lights[index].color != NULL) {
      printf("color: %f, %f, %f\n", lights[index].color[0], lights[index].color[1], lights[index].color[2]);
    }
    if (lights[index].direction != NULL) {
      printf("direction: %f, %f, %f\n", lights[index].direction[0], lights[index].direction[1], lights[index].direction[2]);
    }
    if (lights[index].radialA0 != 0) {
      printf("radialA0: %f\n", lights[index].radialA0);
    }
    if (lights[index].radialA1 != 0) {
      printf("radialA1: %f\n", lights[index].radialA1);
    }
    if (lights[index].radialA2 != 0) {
      printf("radialA2: %f\n", lights[index].radialA2);
    }
    if (lights[index].theta != 0) {
      printf("theta: %f\n", lights[index].theta);
    }
    if (lights[index].angularA0 != 0) {
      printf("angularA0: %f\n", lights[index].angularA0);
    }
    index++;
    printf("\n\n");
  }

}


int main(int argc, char** argv) {

  if (argc != 5){
        printf("usage: raycast width height input.json output.ppm");
        return(1);
    }

    Width = atoi(argv[1]);
    Height = atoi(argv[2]);

    read_scene(argv[3]);

    viewPlane = (Pixel *)malloc(Width * Height * sizeof(Pixel));

    //testPrint();
    raycast();
    write_scene(argv[4], 3);
    return 0;
}
