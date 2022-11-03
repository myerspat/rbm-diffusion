#ifndef _MESH_
#define _MESH_

class Mesh {
private:
  double *_F;
  double **_M;

public:
  Mesh();
  double *getF();
  double **getM();
  ~Mesh();
};

#endif // !_MESH_
