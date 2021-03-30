/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2013 William Oquendo                                   *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

/** @file dem/domain.h .*/

#ifndef MECHSYS_DEM_DOMAIN_H
#define MECHSYS_DEM_DOMAIN_H

// Std lib
#include <cmath>
#include <stdlib.h> // for M_PI
#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include <utility>

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

// Voro++
//#include "src/voro++.cc"
#include "voro++.hh"

// VTK
#ifdef USE_VTK
//#include <vtkCellType.h>
//#include <vtkPolyhedron.h>
#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkPolyData.h>
#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#endif // USE_VTK

// MechSys
#include <mechsys/dem/interacton.h>
#include <mechsys/util/array.h>
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/util/tree.h>

namespace DEM
{

struct MtData;

class Domain
{
public:
    // typedefs
    typedef void (*ptFun_t) (Domain & Dom, void * UserData);

    // Constructor
    Domain(void * UserData=NULL);

    // Destructor
    ~Domain();

    // Particle generation
    void GenSpheres      (int Tag, double L, size_t N, double rho, char const * Type,
                          size_t Randomseed, double fraction, double RminFraction = 1.0);                                        ///< General spheres
    void GenSpheresBox (int Tag, Vec3_t const & X0, Vec3_t const & X1,                                                           ///< Generate spheres within a rectangular box defined by the vectors X0 and X1
                        double R, double rho, char const * Type, size_t Randomseed, double fraction, double RminFraction);
    void GenRice         (int Tag, double L, size_t N, double R, double rho, size_t Randomseed, double fraction);                ///< General rices
    void GenBox          (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf, bool Cohesion=false);            ///< Generate six walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenOpenBox      (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf);                                 ///< Generate five walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenBoundingBox  (int InitialTag, double R, double Cf,bool Cohesion=false);                                              ///< Generate o bounding box enclosing the previous included particles.
    void GenBoundingPlane(int InitialTag, double R, double Cf,bool Cohesion=false);                                              ///< Same as GenBounding but only generates one pair of planes.
    void GenFromMesh     (Mesh::Generic & M, double R, double rho, bool cohesion=false, bool MC=true, double thickness = 0.0);   ///< Generate particles from a FEM mesh generator
    void AddVoroPack     (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    double rho, bool Cohesion, bool Periodic,size_t Randomseed, double fraction, Vec3_t q = OrthoSys::O);                        ///< Generate a Voronoi Packing with dimensions Li and polihedra per side ni
    void AddVoroPack     (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    double rho, bool Cohesion, bVec3_t Periodic,size_t Randomseed, double fraction, Vec3_t q = OrthoSys::O);                     ///< Generate a Voronoi Packing with dimensions Li and polihedra per side ni, Periodic conditions are chosen for each particle
    // Single particle addition
    void AddSphere   (int Tag, Vec3_t const & X, double R, double rho);                                                          ///< Add sphere
    void AddCube     (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddRecBox   (int Tag, Vec3_t const & X, Vec3_t const & L, double R, double rho, double Angle=0, Vec3_t * Axis=NULL);    ///< Add a rectangular box with dimensions given by the vector L
    void AddTetra    (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a tetrahedron at position X with spheroradius R, side of length L and density rho
    void AddDrill    (int Tag, Vec3_t const & X, double R, double Lt, double Ll, double rho);                                    ///< A drill made as a combination of a cube and a pyramid.
    void AddRice     (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a rice at position X with spheroradius R, side of length L and density rho
    void AddPlane    (int Tag, Vec3_t const & X, double R, double Lx,double Ly, double rho, double Angle=0, Vec3_t * Axis=NULL); ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddVoroCell (int Tag, voro::voronoicell & VC, double R, double rho, bool Erode, Vec3_t nv = iVec3_t(1.0,1.0,1.0));      ///< Add a single voronoi cell, it should be built before tough
    void AddTorus    (int Tag, Vec3_t const & X, Vec3_t const & N, double Rmax, double R, double rho);                           ///< Add a single torus at position X with a normal N, circunference Rmax and spheroradius R
    void AddCylinder (int Tag, Vec3_t const & X0, double R0, Vec3_t const & X1, double R1, double R, double rho);                ///< Add a cylinder formed by the connection of two circles at positions X0 and X1 and radii R0 and R1
    void AddFromJson (int Tag, char const * Filename, double R, double rho, double scale,bool Erode = false);                    ///< Add a particle generated from Json mesh


    // 
    //void AddParticle (DEM::Particle * Pa);                                                                                       ///< Add a particle as an exact copy of particle Pa

    // Methods
    void SetProps          (Dict & D);                                                                          ///< Set the properties of individual grains by dictionaries
    void Initialize        (double dt=0.0);                                                                     ///< Set the particles to a initial state and asign the possible insteractions
    void Solve             (double tf, double dt, double dtOut, ptFun_t ptSetup=NULL, ptFun_t ptReport=NULL,
                            char const * FileKey=NULL, size_t VOut=3, size_t Nproc=1,double minEkin=0.0);       ///< Run simulation the simulation up to time tf, with dt and dtOut the time and report steps. The funstion Setup and Report are used to control the workflow form outside, filekey is used to name the report files. VOut has the options 0 no visualization, 1 povray, 2 xmdf and 3 both. minEkin is a minimun of kinetic energy before the simulation stops
    void WritePOV          (char const * FileKey);                                                              ///< Write POV file
    void WriteBPY          (char const * FileKey);                                                              ///< Write BPY (Blender) file
#ifdef USE_HDF5    
    void WriteBF           (char const * FileKey);                                                              ///< Save a h5 with branch and force information
    void WriteFrac         (char const * FileKey);                                                              ///< Save a xdmf file for fracture visualization
    void WriteXDMF         (char const * FileKey);                                                              ///< Save a xdmf file for visualization
    void Save              (char const * FileKey);                                                              ///< Save the current domain
    void Load              (char const * FileKey);                                                              ///< Load the domain form a file
#endif

#ifdef USE_VTK
    void WriteVTKContacts  (char const * FileKey);                                                              ///< Save a vtk - vtp file for conatcs visualization
#endif // USE_VTK

    void UpdateLinkedCells ();                                                                                  ///< Update the linked cells
    void BoundingBox       (Vec3_t & minX, Vec3_t & maxX);                                                      ///< Defines the rectangular box that encloses the particles.
    void Center            (Vec3_t C = Vec3_t(0.0,0.0,0.0));                                                    ///< Centers the domain around C
    void ClearInteractons  ();                                                                                  ///< Reset the interactons
    void ResetInteractons  ();                                                                                  ///< Reset the interactons
    void ResetDisplacements();                                                                                  ///< Reset the displacements
    double MaxDisplacement ();                                                                                  ///< Calculate maximun displacement
    void ResetContacts     ();                                                                                  ///< Reset the displacements
    void ResetBoundaries   ();                                                                                  ///< Reset the Boundary particles
    void EnergyOutput      (size_t IdxOut, std::ostream & OutFile);                                             ///< Output of the energy variables
    void GetGSD            (Array<double> & X, Array<double> & Y, Array<double> & D, size_t NDiv=10) const;     ///< Get the Grain Size Distribution
    void Clusters          ();                                                                                  ///< Check the bounded particles in the domain and how many connected clusters are still present
    void DelParticles      (Array<int> const & Tags);                                                           ///< Delete particle

    // Access methods
    Particle       * GetParticle  (int Tag, bool Check=true);       ///< Find first particle with Tag. Check => check if there are more than one particle with tag=Tag
    Particle const & GetParticle  (int Tag, bool Check=true) const; ///< Find first particle with Tag. Check => check if there are more than one particle with tag=Tag
    void             GetParticles (int Tag, Array<Particle*> & P);  ///< Find all particles with Tag

    // Auxiliar methods
    void   LinearMomentum  (Vec3_t & L);                    ///< Return total momentum of the system
    void   AngularMomentum (Vec3_t & L);                    ///< Return total angular momentum of the system
    double CalcEnergy      (double & Ekin, double & Epot);  ///< Return total energy of the system
    double CriticalDt      ();                              ///< Calculate critical time step from existing particles

#ifdef USE_OMP
    omp_lock_t                                        lck;                         ///< to protect variables in multithreading
#endif
    Array<std::pair<size_t, size_t> >                 ListPosPairs;                ///< List of all possible particles pairs
    iVec3_t                                           LCellDim;                    ///< Dimensions of the linked cell array
    Array<Array <size_t> >                            LinkedCell;                  ///< Linked Cell array for optimization.
    Vec3_t                                            LCxmin;                      ///< Bounding box low   limit for the linked cell array
    Vec3_t                                            LCxmax;                      ///< Bounding box upper limit for the linked cell array

    // Data
    bool                                              Initialized;                 ///< System (particles and interactons) initialized ?
    bool                                              RotPar;                      ///< Check if particles should be rotated, useful if particle displacements are small
    bool                                              Finished;                    ///< Has the simulation finished
    bool                                              Dilate;                      ///< True if eroded particles should be dilated for visualization
    Array<size_t>                                     FreePar;                     ///< Particles that are free
    Array<size_t>                                     NoFreePar;                   ///< Particles that are not free
    Array<Particle*>                                  Particles;                   ///< All particles in domain
    Array<Particle*>                                  ParXmax;                     ///< Particles that are on the Xmax boundary for periodic boudary conditions along the x direction
    Array<Particle*>                                  ParYmax;                     ///< Particles that are on the Ymax boundary for periodic boudary conditions along the y direction
    Array<Particle*>                                  ParXYmax;                    ///< Particles that are on the XYmax boundary for periodic boudary conditions along the y direction
    Array<Interacton*>                                Interactons;                 ///< All interactons
    Array<CInteracton*>                               CInteractons;                ///< Contact interactons
    Array<BInteracton*>                               BInteractons;                ///< Cohesion interactons
    Array<Interacton*>                                PxInteractons;               ///< Interactons for periodic conditions along the x direction
    Array<CInteracton*>                               CPxInteractons;              ///< Contact interacton for periodic conditions along x
    Array<Interacton*>                                PyInteractons;               ///< Interactons for periodic conditions along the y direction
    Array<CInteracton*>                               CPyInteractons;              ///< Contact interacton for periodic conditions along y
    Array<Interacton*>                                PxyInteractons;              ///< Interactons for periodic conditions along the xy direction
    Array<CInteracton*>                               CPxyInteractons;             ///< Contact interacton for periodic conditions along xy
    Vec3_t                                            CamPos;                      ///< Camera position for POV
    double                                            Time;                        ///< Current time
    double                                            Dt;                          ///< Time step
    double                                            Evis;                        ///< Energy dissipated by the viscosity of the grains
    double                                            Efric;                       ///< Energy dissipated by friction
    double                                            Wext;                        ///< Work done by external forces
    double                                            Vs;                          ///< Volume occupied by the grains
    double                                            Ms;                          ///< Total mass of the particles
    double                                            Alpha;                       ///< Verlet distance
    double                                            Beta;                        ///< Binmultiplier
    double                                            Xmax;                        ///< Maximun distance along the X axis (Periodic Boundary)
    double                                            Xmin;                        ///< Minimun distance along the X axis (Periodic Boundary)
    double                                            Ymax;                        ///< Maximun distance along the Y axis (Periodic Boundary)
    double                                            Ymin;                        ///< Minimun distance along the Y axis (Periodic Boundary)
    double                                            MaxDmax;                     ///< Maximun value for the radious of the spheres surronding each particle
    void *                                            UserData;                    ///< Some user data
    String                                            FileKey;                     ///< File Key for output files
    size_t                                            Nproc;                       ///< Number of cores for multithreading
    size_t                                            idx_out;                     ///< Index of output
    std::set<std::pair<Particle *, Particle *> >      Listofpairs;                 ///< List of pair of particles associated per interacton for memory optimization
    std::set<std::pair<Particle *, Particle *> >      PxListofpairs;               ///< List of pair of particles associated per interacton for memory optimization under periodic boundary conditions
    std::set<std::pair<Particle *, Particle *> >      PyListofpairs;               ///< List of pair of particles associated per interacton for memory optimization under periodic boundary conditions
    std::set<std::pair<Particle *, Particle *> >      PxyListofpairs;              ///< List of pair of particles associated per interacton for memory optimization under periodic boundary conditions
    Array<Array <int> >                               Listofclusters;              ///< List of particles belonging to bounded clusters (applies only for cohesion simulations)
    MtData *                                          MTD;                         ///< Multithread data

    // Some utilities when the interactions are mainly between spheres
    bool                                              MostlySpheres;               ///< If the simulation is mainly between spheres this should be true
    FrictionMap_t                                     FricSpheres;                 ///< The friction value for spheres only
    FrictionMap_t                                     RollSpheres;                 ///< Map storing the rolling resistance between spheres
    void     CalcForceSphere();                                                    ///< Calculate force between only spheres spheres
    
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


struct MtData   /// A structure for the multi-thread data
{
    size_t                       ProcRank; ///< Rank of the thread
    size_t                         N_Proc; ///< Total number of threads
    DEM::Domain *                     Dom; ///< Pointer to the lbm domain
    double                            Dmx; ///< Maximun displacement
    Array<std::pair<size_t,size_t> >   LC; ///< A temporal list of new contacts
    Array<size_t>                     LCI; ///< A temporal array of posible Cinteractions
    Array<size_t>                     LCB; ///< A temporal array of posible Binteractions
    Array<std::pair<size_t,size_t> >  LPC; ///< A temporal list of new contacts for periodic boundary conditions
    Array<size_t>                    LPCI; ///< A temporal array of posible Cinteractions for periodic boundary conditions
    Array<size_t>                     LBP; ///< A temporal array of possible boundary particles
    Array<std::pair<iVec3_t,size_t> > LLC; ///< A temporal array of possible linked cells locations
    Array<std::pair<size_t,size_t> >  LPP; ///< A temporal array of possible partcle types
};

// Constructor & Destructor

inline Domain::Domain (void * UD)
    :  Initialized(false), Dilate(false), RotPar(true), Time(0.0), Alpha(0.05), Beta(1.0), UserData(UD)
{
    MostlySpheres = false;
    Xmax = Xmin = Ymax = Ymin = 0.0;
    CamPos = 1.0, 2.0, 3.0;
#ifdef USE_OMP
    omp_init_lock(&lck);
#endif
}

inline Domain::~Domain ()
{
    for (size_t i=0; i<Particles.Size();   ++i) if (Particles  [i]!=NULL) delete Particles  [i];
    for (size_t i=0; i<Interactons.Size(); ++i) if (Interactons[i]!=NULL) delete Interactons[i];
}

//All the methods for particle generation

#include<mechsys/dem/dompargen.h>

// Methods

inline void Domain::SetProps (Dict & D)
{
    for (size_t i =0 ; i<Particles.Size(); i++)
    {
        for (size_t j=0; j<D.Keys.Size(); ++j)
        {
            int tag = D.Keys[j];
            if (tag==Particles[i]->Tag)
            {
                SDPair const & p = D(tag);
                if (p.HasKey("Gn"))
                {
                    Particles[i]->Props.Gn = p("Gn");
                }
                if (p.HasKey("Gt"))
                {
                    Particles[i]->Props.Gt = p("Gt");
                }
                if (p.HasKey("Gv"))
                {
                    Particles[i]->Props.Gv = p("Gv");
                }
                if (p.HasKey("Gm"))
                {
                    Particles[i]->Props.Gm = p("Gm");
                }
                if (p.HasKey("Kn"))
                {
                    Particles[i]->Props.Kn = p("Kn");
                }
                if (p.HasKey("Kt"))
                {
                    Particles[i]->Props.Kt = p("Kt");
                }
                if (p.HasKey("Bn"))
                {
                    Particles[i]->Props.Bn = p("Bn");
                }
                if (p.HasKey("Bt"))
                {
                    Particles[i]->Props.Bt = p("Bt");
                }
                if (p.HasKey("Bm"))
                {
                    Particles[i]->Props.Bm = p("Bm");
                }
                if (p.HasKey("Mu"))
                {
                    Particles[i]->Props.Mu = p("Mu");
                }
                if (p.HasKey("Eps"))
                {
                    Particles[i]->Props.eps = p("Eps");
                }
                if (p.HasKey("Beta"))
                {
                    Particles[i]->Props.Beta = p("Beta");
                }
                if (p.HasKey("Eta"))
                {
                    Particles[i]->Props.Eta = p("Eta");
                }
            }
        }
    }
    for (size_t i=0; i<BInteractons.Size(); i++)
    {
        BInteractons[i]->UpdateParameters();
    }
    for (size_t i=0; i<CInteractons.Size(); i++)
    {
        CInteractons[i]->UpdateParameters();
    }
}

inline void Domain::Initialize (double dt)
{
    if (!Initialized)
    {
        // initialize all particles
        for (size_t i=0; i<Particles.Size(); i++)
        {
            Particles[i]->Initialize(i);
            Particles[i]->InitializeVelocity(dt);
        }
        //Initializing the energies
        Evis = 0.0;
        Efric = 0.0;
        Wext = 0.0;

        // initialize
        //ResetInteractons();
        // info
        Util::Stopwatch stopwatch;
        printf("\n%s--- Initializing particles ------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
        // set flag
        Initialized = true;

        // info
        double Ekin, Epot, Etot;
        Etot = CalcEnergy (Ekin, Epot);
        printf("%s  Kinematic energy   = %g%s\n",TERM_CLR4, Ekin, TERM_RST);
        printf("%s  Potential energy   = %g%s\n",TERM_CLR4, Epot, TERM_RST);
        printf("%s  Total energy       = %g%s\n",TERM_CLR2, Etot, TERM_RST);
    }
    else
    {
        for (size_t i=0; i<Particles.Size(); i++)
        {
            if (Particles[i]->vxf) Particles[i]->xb(0) = Particles[i]->x(0) - Particles[i]->v(0)*dt;
            if (Particles[i]->vyf) Particles[i]->xb(1) = Particles[i]->x(1) - Particles[i]->v(1)*dt;
            if (Particles[i]->vzf) Particles[i]->xb(2) = Particles[i]->x(2) - Particles[i]->v(2)*dt;
        }
    }

}

inline void Domain::Solve (double tf, double dt, double dtOut, ptFun_t ptSetup, ptFun_t ptReport, char const * TheFileKey, size_t VOut, size_t TheNproc, double minEkin)
{
    if (VOut > 3) throw new Fatal("Domain::Solve The visualization argument can only have 4 values: 0 None, 1 povray visualization, 2 xdmf visualization and 3 both options");
    // Assigning some domain particles especifically to the output
    FileKey.Printf("%s",TheFileKey);
    idx_out = 0;
    
    //Assigning the vlaue for the time step to the domain variable
    Dt = dt;
    Nproc = TheNproc;

    // initialize particles
    Initialize (Dt);


    // calc the total volume of particles (solids)
    FreePar.Resize(0);
    NoFreePar.Resize(0);
    Vs = 0.0;
    Ms = 0.0;
    MaxDmax        =  0.0;
    double MaxKn   =  0.0;
    double MaxBn   =  0.0;
    double MinDmax = -1.0;
    double MinMass = -1.0;
    for (size_t i=0; i<Particles.Size(); i++) 
    { 
        if (Particles[i]->IsFree())
        {
            Vs += Particles[i]->Props.V;
            Ms += Particles[i]->Props.m;
            if (Particles[i]->Dmax     > MaxDmax) MaxDmax = Particles[i]->Dmax;
            if (Particles[i]->Props.Kn > MaxKn  ) MaxKn   = Particles[i]->Props.Kn;
            if (Particles[i]->Dmax     < MinDmax||(MinDmax<0.0)) MinDmax = Particles[i]->Dmax;
            if (Particles[i]->Props.m  < MinMass||(MinMass<0.0)) MinMass = Particles[i]->Props.m;
            FreePar.Push(i);
        }
        else NoFreePar.Push(i);
    }
    for (size_t i=0; i<BInteractons.Size(); i++)
    {
        double pbn = std::max(BInteractons[i]->Bn/BInteractons[i]->L0,BInteractons[i]->Bt/BInteractons[i]->L0);
        if (pbn > MaxBn) MaxBn = pbn;
    }


    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    printf("%s  Total mass   of free particles            =  %g%s\n"        ,TERM_CLR4, Ms                                   , TERM_RST);
    printf("%s  Total volume of free particles            =  %g%s\n"        ,TERM_CLR4, Vs                                   , TERM_RST);
    printf("%s  Total number of particles                 =  %zd%s\n"       ,TERM_CLR2, Particles.Size()                     , TERM_RST);
    printf("%s  Time step                                 =  %g%s\n"        ,TERM_CLR2, dt                                   , TERM_RST);
    printf("%s  Verlet distance                           =  %g%s\n"        ,TERM_CLR2, Alpha                                , TERM_RST);
    printf("%s  Suggested Time Step                       =  %g%s\n"        ,TERM_CLR5, 0.1*sqrt(MinMass/(MaxKn+MaxBn))      , TERM_RST);
    printf("%s  Suggested Verlet distance                 =  %g or  %g%s\n" ,TERM_CLR5, 0.5*MinDmax, 0.25*(MinDmax + MaxDmax), TERM_RST);
    if (fabs(Xmax-Xmin)>1.0e-12)
    printf("%s  Periodic Boundary conditions in X between =  %g and %g%s\n" ,TERM_CLR5, Xmin, Xmax                           , TERM_RST);
    if (fabs(Ymax-Ymin)>1.0e-12)
    printf("%s  Periodic Boundary conditions in Y between =  %g and %g%s\n" ,TERM_CLR5, Ymin, Ymax                           , TERM_RST);

    if (Alpha > Beta*MaxDmax)
    {
        Alpha = Beta*MaxDmax;
        printf("%s  Verlet distance changed to       =  %g%s\n"   ,TERM_CLR2, Alpha                                    , TERM_RST);
    }

    // solve
    double t0   = Time;     // initial time
    double tout = t0; // time position for output

    Finished = false;

    // string to output energy data, if user gives the FileKey
    std::ostringstream oss_energy; 
    EnergyOutput (idx_out, oss_energy);

    MTD = new DEM::MtData[Nproc];
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].N_Proc   = Nproc;
        MTD[i].ProcRank = i;
        MTD[i].Dom      = this;
        MTD[i].Dmx      = 0.0;
    }
    //for (size_t i=0; i<Particles.Size()-1; i++)
    //for (size_t j=i+1; j<Particles.Size(); j++)
    //{
        //ListPosPairs.Push(std::make_pair(i,j));
    //}
    LinkedCell.Resize(0);
    BoundingBox(LCxmin,LCxmax);
    LCellDim = (LCxmax - LCxmin)/(2.0*Beta*MaxDmax) + iVec3_t(1,1,1);
    //LCellDim = (LCxmax - LCxmin)/(2.0*Beta*MaxDmax);
    LinkedCell.Resize(LCellDim(0)*LCellDim(1)*LCellDim(2));
    //std::cout << LCellDim << " " << MaxDmax << std::endl;
    //iVec3_t iv;
    //idx2Pt(16,iv,LCellDim);
    //std::cout << iv << std::endl;
#ifdef USE_OMP
    //std::cout << "2 " << CInteractons.Size() << std::endl;
    //ResetDisplacements
    ResetDisplacements();

    //std::cout << "3 " << CInteractons.Size() << std::endl;
    //UpdateLinkedCells
    UpdateLinkedCells();

    //std::cout << ListPosPairs.Size() << endl;
    //std::cout << "4 " << CInteractons.Size() << std::endl;
    //ResetContacts
    ResetContacts();
#else

    // set the displacement of the particles to zero (for the Halo)
    ResetDisplacements();

    // build the map of possible contacts (for the Halo)
    ResetContacts();
    if (fabs(Xmax-Xmin)>Alpha) ResetBoundaries();
    if (fabs(Ymax-Ymin)>Alpha) ResetBoundaries();

#endif

    // run
    while (Time<tf)
    {

        // output
        if (Time>=tout)
        {
            double Ekin,Epot;
            CalcEnergy(Ekin,Epot);
            if (Ekin<minEkin&&Time>0.1*tf)
            {
                printf("\n%s--- Minimun energy reached ---------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
                break;
            }
            if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            if (TheFileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
#ifdef USE_HDF5
                if (VOut==2||VOut==3)    WriteXDMF    (fn.CStr());
#endif
                if (VOut==1||VOut==3)    WritePOV     (fn.CStr());
                //EnergyOutput (idx_out, oss_energy);
            }
            if (BInteractons.Size()>0) Clusters();
            idx_out++;
            tout += dtOut;
        }
#ifdef USE_OMP 
        //Initialize particles
        //std::cout << "1" << std::endl;
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0; i<Particles.Size(); i++)
        {
            // set the force and torque to the fixed values
            Particles[i]->F = Particles[i]->Ff;
            Particles[i]->T = Particles[i]->Tf;

            //Particles[i]->Bdry = false;
        }

        //Calculate forces
        //std::cout << Interactons.Size() << " 2" << std::endl;
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0; i<Interactons.Size(); i++)
        {
		    if (Interactons[i]->CalcForce(Dt))
            {
                String f_error(FileKey+"_error");
                Save     (f_error.CStr());
                WriteXDMF(f_error.CStr());
                std::cout << "Maximun overlap detected between particles at time " << Time << std::endl;
                sleep(1);
                throw new Fatal("Maximun overlap detected between particles");
            }
            omp_set_lock  (&Interactons[i]->P1->lck);
            Interactons[i]->P1->F += Interactons[i]->F1;
            Interactons[i]->P1->T += Interactons[i]->T1;
            omp_unset_lock(&Interactons[i]->P1->lck);
            omp_set_lock  (&Interactons[i]->P2->lck);
            Interactons[i]->P2->F += Interactons[i]->F2;
            Interactons[i]->P2->T += Interactons[i]->T2;
            omp_unset_lock(&Interactons[i]->P2->lck);
        }

        if(MostlySpheres) CalcForceSphere();
        // Periodic Boundary
        //std::cout << "3" << std::endl;
        if (Xmax-Xmin>Alpha)
        {
            Vec3_t v(Xmin-Xmax,0.0,0.0);
            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<ParXmax.Size(); i++) ParXmax[i]->Translate(v);

            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<PxInteractons.Size(); i++)
            {
		        if (PxInteractons[i]->CalcForce(Dt))
                {
                    String f_error(FileKey+"_error");
                    Save     (f_error.CStr());
                    WriteXDMF(f_error.CStr());
                    std::cout << "Maximun overlap detected between particles at the periodic X boundary at time " << Time << std::endl;
                    std::cout << ParYmax.Has(PxInteractons[i]->P1) << std::endl;
                    std::cout << ParYmax.Has(PxInteractons[i]->P2) << std::endl;
                    std::cout << ParXmax.Has(PxInteractons[i]->P1) << std::endl;
                    std::cout << ParXmax.Has(PxInteractons[i]->P2) << std::endl;
                    sleep(1);
                    throw new Fatal("Maximun overlap detected between particles");
                }
                omp_set_lock  (&PxInteractons[i]->P1->lck);
                PxInteractons[i]->P1->F += PxInteractons[i]->F1;
                PxInteractons[i]->P1->T += PxInteractons[i]->T1;
                omp_unset_lock(&PxInteractons[i]->P1->lck);
                omp_set_lock  (&PxInteractons[i]->P2->lck);
                PxInteractons[i]->P2->F += PxInteractons[i]->F2;
                PxInteractons[i]->P2->T += PxInteractons[i]->T2;
                omp_unset_lock(&PxInteractons[i]->P2->lck);
            }

            v = Vec3_t(Xmax-Xmin,0.0,0.0);
            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<ParXmax.Size(); i++) ParXmax[i]->Translate(v);
        }
        
        if (Ymax-Ymin>Alpha)
        {
            Vec3_t v(0.0,Ymin-Ymax,0.0);
            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<ParYmax.Size(); i++) ParYmax[i]->Translate(v);

            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<PyInteractons.Size(); i++)
            {
		        if (PyInteractons[i]->CalcForce(Dt))
                {
                    String f_error(FileKey+"_error");
                    Save     (f_error.CStr());
                    WriteXDMF(f_error.CStr());
                    std::cout << "Maximun overlap detected between particles at the periodic Y boundary at time " << Time << std::endl;
                    sleep(1);
                    throw new Fatal("Maximun overlap detected between particles");
                }
                omp_set_lock  (&PyInteractons[i]->P1->lck);
                PyInteractons[i]->P1->F += PyInteractons[i]->F1;
                PyInteractons[i]->P1->T += PyInteractons[i]->T1;
                omp_unset_lock(&PyInteractons[i]->P1->lck);
                omp_set_lock  (&PyInteractons[i]->P2->lck);
                PyInteractons[i]->P2->F += PyInteractons[i]->F2;
                PyInteractons[i]->P2->T += PyInteractons[i]->T2;
                omp_unset_lock(&PyInteractons[i]->P2->lck);
            }

            v = Vec3_t(0.0,Ymax-Ymin,0.0);
            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<ParYmax.Size(); i++) ParYmax[i]->Translate(v);
        }
        
        if ((Ymax-Ymin>Alpha)&&(Xmax-Xmin>Alpha))
        {
            Vec3_t v(Xmin-Xmax,Ymin-Ymax,0.0);
            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<ParXYmax.Size(); i++) ParXYmax[i]->Translate(v);

            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<PxyInteractons.Size(); i++)
            {
		        if (PxyInteractons[i]->CalcForce(Dt))
                {
                    String f_error(FileKey+"_error");
                    Save     (f_error.CStr());
                    WriteXDMF(f_error.CStr());
                    std::cout << "Maximun overlap detected between particles at the periodic XY boundary at time " << Time << std::endl;
                    sleep(1);
                    throw new Fatal("Maximun overlap detected between particles");
                }
                omp_set_lock  (&PxyInteractons[i]->P1->lck);
                PxyInteractons[i]->P1->F += PxyInteractons[i]->F1;
                PxyInteractons[i]->P1->T += PxyInteractons[i]->T1;
                omp_unset_lock(&PxyInteractons[i]->P1->lck);
                omp_set_lock  (&PxyInteractons[i]->P2->lck);
                PxyInteractons[i]->P2->F += PxyInteractons[i]->F2;
                PxyInteractons[i]->P2->T += PxyInteractons[i]->T2;
                omp_unset_lock(&PxyInteractons[i]->P2->lck);
            }

            v = Vec3_t(Xmax-Xmin,Ymax-Ymin,0.0);
            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<ParXYmax.Size(); i++) ParXYmax[i]->Translate(v);
        }
        // tell the user function to update its data
        //std::cout << "4" << std::endl;
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);

        // Move Particles
        //std::cout << "5" << std::endl;
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<Nproc;i++)
        {
            MTD[i].Dmx = 0.0;
        }

        if (RotPar)
        {
            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<Particles.Size(); i++)
            {
                //std::cout << "1" << std::endl;
		        Particles[i]->Translate(Dt);
                //std::cout << "2" << std::endl;
		        Particles[i]->Rotate(Dt);
                //std::cout << "3" << std::endl;
                if (Particles[i]->MaxDisplacement()>MTD[omp_get_thread_num()].Dmx) MTD[omp_get_thread_num()].Dmx = Particles[i]->MaxDisplacement();
            }
        }
        else
        {
            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<Particles.Size(); i++)
            {
                //std::cout << "1" << std::endl;
		        Particles[i]->Translate(Dt);
                //std::cout << "2" << std::endl;
                if (Particles[i]->MaxDisplacement()>MTD[omp_get_thread_num()].Dmx) MTD[omp_get_thread_num()].Dmx = Particles[i]->MaxDisplacement();
            }
        }

        double maxdis = 0.0;
        for (size_t i=0;i<Nproc;i++)
        {
            if (maxdis<MTD[i].Dmx) maxdis = MTD[i].Dmx;
        }

        //Update Linked Cells
        //std::cout << "6" << std::endl;
        if (maxdis>Alpha)
        {
            //std::cout << "A" <<  std::endl;
            LinkedCell.Resize(0);
            BoundingBox(LCxmin,LCxmax);
            LCellDim = (LCxmax - LCxmin)/(2.0*Beta*MaxDmax) + iVec3_t(1,1,1);
            //LCellDim = (LCxmax - LCxmin)/(2.0*Beta*MaxDmax);
            LinkedCell.Resize(LCellDim(0)*LCellDim(1)*LCellDim(2));

            //ResetDisplacements
            ResetDisplacements();

            //UpdateLinkedCells
            UpdateLinkedCells();

            //ResetContacts
            ResetContacts();
        }

#else 
       //Not used anymore, OpenMP should be used by default.

#endif
        

        // next time position
        Time += Dt;
    }

    // last output
    Finished = true;
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);

    // save energy data
    //if (TheFileKey!=NULL)
    //{
        //String fn;
        //fn.Printf("%s_energy.res",TheFileKey);
        //std::ofstream fe(fn.CStr());
        //fe << oss_energy.str();
        //fe.close();
    //}

    // info
    double Ekin, Epot, Etot;
    Etot = CalcEnergy (Ekin, Epot);
    printf("%s  Kinematic energy   = %g%s\n",TERM_CLR4, Ekin, TERM_RST);
    printf("%s  Potential energy   = %g%s\n",TERM_CLR4, Epot, TERM_RST);
    printf("%s  Total energy       = %g%s\n",TERM_CLR2, Etot, TERM_RST);
}

inline void Domain::WritePOV (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".pov");
    std::ofstream of(fn.CStr(), std::ios::out);
    POVHeader (of);
    POVSetCam (of, CamPos, OrthoSys::O);
    Array <String> Colors(10);
    Colors = "Gray","Blue","Yellow","Gold","Green","Blue","Orange","Salmon","Copper","Aquamarine";
    for (size_t i=0; i<Particles.Size(); i++)
    {
        if (!Particles[i]->IsFree()) Particles[i]->Draw(of,"Col_Glass_Bluish");
        else
        {
            bool found = false;
            for (size_t j=0;j<Listofclusters.Size();j++)
            {
                if (Listofclusters[j].Has(i)&&BInteractons.Size()>0)
                {
                    Particles[i]->Draw(of,Colors[j%10].CStr());
                    found = true;
                    break;
                }
            }
            if (!found) Particles[i]->Draw(of,"Red");
        }
    }

}

inline void Domain::WriteBPY (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".bpy");
    std::ofstream of(fn.CStr(), std::ios::out);
    BPYHeader(of);
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Draw (of,"",true);
    of.close();
}

#ifdef USE_HDF5

inline void Domain::WriteBF (char const * FileKey)
{

    size_t n_fn = 0;
    size_t n_rl = 0;

    for (size_t i=0;i<CInteractons.Size();i++)
    {
        //if ((norm(CInteractons[i]->Fnet)>1.0e-12)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree())) n_fn++;
        if (norm(CInteractons[i]->Fnet)>1.0e-12)
        {
            n_fn++;
            if (CInteractons[i]->P1->Verts.Size()==1&&CInteractons[i]->P2->Verts.Size()==1) n_rl++;
        }
    }

    if (n_fn==0) return;
    
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    float  *  Fnnet = new float[3*n_fn];
    float  *  Ftnet = new float[3*n_fn];
    float  *  Froll = new float[3*n_fn];
    float  * Branch = new float[3*n_fn];
    float  *   Orig = new float[3*n_fn];
    int    *    ID1 = new   int[  n_fn];
    int    *    ID2 = new   int[  n_fn];

    size_t idx = 0;

    // Saving Collision forces
    for (size_t i=0;i<CInteractons.Size();i++)
    {
        //if ((norm(CInteractons[i]->Fnet)>1.0e-12)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree()))
        if (norm(CInteractons[i]->Fnet)>1.0e-12)
        {
            Fnnet [3*idx  ] = float(CInteractons[i]->Fnet  (0));
            Fnnet [3*idx+1] = float(CInteractons[i]->Fnet  (1));
            Fnnet [3*idx+2] = float(CInteractons[i]->Fnet  (2));
            Ftnet [3*idx  ] = float(CInteractons[i]->Ftnet (0));
            Ftnet [3*idx+1] = float(CInteractons[i]->Ftnet (1));
            Ftnet [3*idx+2] = float(CInteractons[i]->Ftnet (2));
            if (n_rl>0)
            {
            Froll [3*idx  ] = float(CInteractons[i]->Fn    (0));
            Froll [3*idx+1] = float(CInteractons[i]->Fn    (1));
            Froll [3*idx+2] = float(CInteractons[i]->Fn    (2));
            }
            Branch[3*idx  ] = float(CInteractons[i]->P1->x(0)-CInteractons[i]->P2->x(0));
            Branch[3*idx+1] = float(CInteractons[i]->P1->x(1)-CInteractons[i]->P2->x(1)); 
            Branch[3*idx+2] = float(CInteractons[i]->P1->x(2)-CInteractons[i]->P2->x(2)); 
            //Orig  [3*idx  ] = 0.5*float(CInteractons[i]->P1->x(0)+CInteractons[i]->P2->x(0));
            //Orig  [3*idx+1] = 0.5*float(CInteractons[i]->P1->x(1)+CInteractons[i]->P2->x(1)); 
            //Orig  [3*idx+2] = 0.5*float(CInteractons[i]->P1->x(2)+CInteractons[i]->P2->x(2)); 
            Orig  [3*idx  ] = float(CInteractons[i]->P2->x(0));
            Orig  [3*idx+1] = float(CInteractons[i]->P2->x(1)); 
            Orig  [3*idx+2] = float(CInteractons[i]->P2->x(2)); 
            ID1   [idx]     = int  (CInteractons[i]->P1->Index);
            ID2   [idx]     = int  (CInteractons[i]->P2->Index);
            idx++;
        }
    }

    hsize_t dims[1];
    dims[0] = 3*n_fn;
    String dsname;
    dsname.Printf("Normal");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Fnnet );
    dsname.Printf("Tangential");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ftnet );
    if (n_rl>0)
    {
    dsname.Printf("Rolling");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Froll );
    }
    dsname.Printf("Branch");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Branch);
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Orig);
    dims[0] = n_fn;
    dsname.Printf("ID1");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,ID1   );
    dsname.Printf("ID2");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,ID2   );


    delete [] Fnnet;
    delete [] Ftnet;
    delete [] Froll;
    delete [] Branch;
    delete [] Orig;
    delete [] ID1;
    delete [] ID2;

    //Saving Cohesive forces
    if (BInteractons.Size()>0)
    {
    float  *  Bnnet = new float[3*BInteractons.Size()];
    float  *  Btnet = new float[3*BInteractons.Size()];
    float  *  BOrig = new float[3*BInteractons.Size()];
    int    *  BVal  = new   int[  BInteractons.Size()];
    int    *  BID1  = new   int[  BInteractons.Size()];
    int    *  BID2  = new   int[  BInteractons.Size()];

    idx = 0;
    for (size_t i=0;i<BInteractons.Size();i++)
    {
        Bnnet [3*idx  ] = float(BInteractons[i]->Fnet  (0));
        Bnnet [3*idx+1] = float(BInteractons[i]->Fnet  (1));
        Bnnet [3*idx+2] = float(BInteractons[i]->Fnet  (2));
        Btnet [3*idx  ] = float(BInteractons[i]->Ftnet (0));
        Btnet [3*idx+1] = float(BInteractons[i]->Ftnet (1));
        Btnet [3*idx+2] = float(BInteractons[i]->Ftnet (2));
        BOrig [3*idx  ] = float(BInteractons[i]->xnet(0));
        BOrig [3*idx+1] = float(BInteractons[i]->xnet(1)); 
        BOrig [3*idx+2] = float(BInteractons[i]->xnet(2)); 
        BVal  [idx]     = int  (BInteractons[i]->valid);
        BID1  [idx]     = int  (BInteractons[i]->P1->Index);
        BID2  [idx]     = int  (BInteractons[i]->P2->Index);
        idx++;

    }
    hsize_t dims[1];
    dims[0] = 3*BInteractons.Size();
    String dsname;
    dsname.Printf("BNormal");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Bnnet );
    dsname.Printf("BTangential");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Btnet );
    dsname.Printf("BPosition");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,BOrig );
    dims[0] = BInteractons.Size();
    dsname.Printf("BVal");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,BVal  );
    dsname.Printf("BID1");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,BID1  );
    dsname.Printf("BID2");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,BID2  );

    delete [] Bnnet;
    delete [] Btnet;
    delete [] BOrig;
    delete [] BVal ;
    delete [] BID1 ;
    delete [] BID2 ;


    }


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"BranchForce\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << n_fn << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << n_fn << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Normal\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Normal \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tangential\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tangential \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    if (n_rl>0)
    {
    oss << "     <Attribute Name=\"Rolling\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Rolling \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    }
    oss << "     <Attribute Name=\"Branch\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Branch \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"ID1\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/ID1 \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"ID2\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/ID2 \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    if (BInteractons.Size()>0)
    {
    oss << "   <Grid Name=\"CohesiveForce\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << BInteractons.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << BInteractons.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/BPosition \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"BNormal\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BNormal \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"BTangential\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BTangential \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"BVal\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BVal \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"BID1\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BID1 \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"BID2\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BID2 \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    }
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Domain::WriteXDMF (char const * FileKey)
{
    size_t N_Faces = 0;
    size_t N_Verts = 0;
    for (size_t i=0; i<Particles.Size(); i++) 
    { 
        for (size_t j=0;j<Particles[i]->Faces.Size();j++)
        {
            N_Faces += Particles[i]->Faces[j]->Edges.Size();
        }
        N_Verts += Particles[i]->Verts.Size() + Particles[i]->Faces.Size();
    }

    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if (N_Faces>0)
    {

        //Geometric information
        float  * Verts   = new float [3*N_Verts];
        int    * FaceCon = new int   [3*N_Faces];
        
        //Atributes
        int    * Tags    = new int   [  N_Faces];
        int    * Clus    = new int   [  N_Faces];
        float  * Vel     = new float [  N_Faces];
        float  * Ome     = new float [  N_Faces];
        //float  * Stress  = new float [9*N_Faces];

        size_t n_verts = 0;
        size_t n_faces = 0;
        size_t n_attrs = 0;
        //size_t n_attrv = 0;
        //size_t n_attrt = 0;
        for (size_t i=0;i<Particles.Size();i++)
        {
            Particle * Pa = Particles[i];
            size_t n_refv = n_verts/3;
            Array<Vec3_t> Vtemp(Pa->Verts.Size());
            Array<Vec3_t> Vres (Pa->Verts.Size());
            for (size_t j=0;j<Pa->Verts.Size();j++)
            {
                Vtemp[j] = *Pa->Verts[j];
                Vres [j] = *Pa->Verts[j];
            }
            double multiplier = 0.0;
            if (Dilate&&Pa->Eroded&&Pa->Faces.Size()>=4)
            {
                DEM::Dilation(Vtemp,Pa->EdgeCon,Pa->FaceCon,Vres,Pa->Props.R);
                multiplier = 1.0;
            }
            for (size_t j=0;j<Pa->Verts.Size();j++)
            {
                //Verts[n_verts++] = (float) (*Pa->Verts[j])(0);
                //Verts[n_verts++] = (float) (*Pa->Verts[j])(1);
                //Verts[n_verts++] = (float) (*Pa->Verts[j])(2);
                Verts[n_verts++] = float(Vres[j](0));
                Verts[n_verts++] = float(Vres[j](1));
                Verts[n_verts++] = float(Vres[j](2));
            }
            size_t n_reff = n_verts/3;
            for (size_t j=0;j<Pa->FaceCon.Size();j++)
            {
                Vec3_t C,N;
                Pa->Faces[j]->Centroid(C);
                Pa->Faces[j]->Normal(N);
                Verts[n_verts++] = float(C(0) + multiplier*Pa->Props.R*N(0));
                Verts[n_verts++] = float(C(1) + multiplier*Pa->Props.R*N(1));
                Verts[n_verts++] = float(C(2) + multiplier*Pa->Props.R*N(2));
                //Verts[n_verts++] = (float) C(0);
                //Verts[n_verts++] = (float) C(1);
                //Verts[n_verts++] = (float) C(2);
                for (size_t k=0;k<Pa->FaceCon[j].Size();k++)
                {
                    size_t nin = Pa->FaceCon[j][k];
                    size_t nen = Pa->FaceCon[j][(k+1)%Pa->FaceCon[j].Size()];
                    FaceCon[n_faces++] = int(n_reff + j);  
                    FaceCon[n_faces++] = int(n_refv + nin);
                    FaceCon[n_faces++] = int(n_refv + nen);

                    //Writing the attributes
                    Tags  [n_attrs] = int(Pa->Tag);
                    Clus  [n_attrs] = size_t(Pa->Cluster);
                    Vel   [n_attrs] = float(norm(Pa->v));
                    Ome   [n_attrs] = float(norm(Pa->w));
                    n_attrs++;

                    //Vel [n_attrv  ] = (float) Pa->v(0);
                    //Vel [n_attrv+1] = (float) Pa->v(1);
                    //Vel [n_attrv+2] = (float) Pa->v(2);
                    //n_attrv += 3;

                    //Stress[n_attrt  ] = (float) Pa->M(0,0);
                    //Stress[n_attrt+1] = (float) Pa->M(1,0);
                    //Stress[n_attrt+2] = (float) Pa->M(2,0);
                    //Stress[n_attrt+3] = (float) Pa->M(0,1);
                    //Stress[n_attrt+4] = (float) Pa->M(1,1);
                    //Stress[n_attrt+5] = (float) Pa->M(2,1);
                    //Stress[n_attrt+6] = (float) Pa->M(0,2);
                    //Stress[n_attrt+7] = (float) Pa->M(1,2);
                    //Stress[n_attrt+8] = (float) Pa->M(2,2);
                    //n_attrt += 9;
                }
            }
        }

        //Write the data
        hsize_t dims[1];
        String dsname;
        dims[0] = 3*N_Verts;
        dsname.Printf("Verts");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Verts);
        dims[0] = 3*N_Faces;
        dsname.Printf("FaceCon");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,FaceCon);
        dims[0] = N_Faces;
        dsname.Printf("Tag");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Tags   );
        dims[0] = N_Faces;
        dsname.Printf("Cluster");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Clus   );
        dims[0] = N_Faces;
        dsname.Printf("Velocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Vel);
        dims[0] = N_Faces;
        dsname.Printf("AngVelocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ome);
        
        //Erasing the data
        delete [] Verts;
        delete [] FaceCon;
        delete [] Tags;
        delete [] Clus;
        delete [] Vel;
        delete [] Ome;
        //delete [] Stress;
    }
    // Storing center of mass data
    
    float * Radius = new float[  Particles.Size()];
    float * Posvec = new float[3*Particles.Size()];
    float * Velvec = new float[3*Particles.Size()];
    float * Omevec = new float[3*Particles.Size()];
    float * Amovec = new float[3*Particles.Size()];
    float * Ekin   = new float[  Particles.Size()];
    int   * Tag    = new int  [  Particles.Size()];

    for (size_t i=0;i<Particles.Size();i++)
    {
        Vec3_t Ome,L,t1,t2;
        Rotation(Particles[i]->w,Particles[i]->Q,Ome);
        t1 = Particles[i]->I(0)*Particles[i]->w(0),Particles[i]->I(1)*Particles[i]->w(1),Particles[i]->I(2)*Particles[i]->w(2);
        Rotation (t1,Particles[i]->Q,t2);
        L = Particles[i]->Props.m*cross(Particles[i]->x,Particles[i]->v)+t2;


        Particles[i]->Verts.Size()==1 ? Radius[i] = float(Particles[i]->Dmax) : Radius[i] = 0.0;
        Posvec[3*i  ] = float(Particles[i]->x(0));
        Posvec[3*i+1] = float(Particles[i]->x(1));
        Posvec[3*i+2] = float(Particles[i]->x(2));
        Velvec[3*i  ] = float(Particles[i]->v(0));
        Velvec[3*i+1] = float(Particles[i]->v(1));
        Velvec[3*i+2] = float(Particles[i]->v(2));
        Omevec[3*i  ] = float(Ome(0));
        Omevec[3*i+1] = float(Ome(1)); 
        Omevec[3*i+2] = float(Ome(2)); 
        Amovec[3*i  ] = float(L(0));
        Amovec[3*i+1] = float(L(1)); 
        Amovec[3*i+2] = float(L(2)); 
        Ekin  [i]     = float(Particles[i]->Ekin+Particles[i]->Erot);
        Tag   [i]     = int  (Particles[i]->Tag);  
    }

    hsize_t dims[1];
    dims[0] = 3*Particles.Size();
    String dsname;
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    dsname.Printf("PVelocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    dsname.Printf("PAngVelocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Omevec);
    dsname.Printf("PAngMomentum");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Amovec);
    dims[0] = Particles.Size();
    dsname.Printf("Radius");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Radius);
    dsname.Printf("PEkin");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ekin);
    dsname.Printf("PTag");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,Tag   );


    delete [] Radius;
    delete [] Posvec;
    delete [] Velvec;
    delete [] Omevec;
    delete [] Amovec;
    delete [] Ekin;
    delete [] Tag;


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);
    

    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    if(N_Faces>0)
    {
    oss << "   <Grid Name=\"DEM_Faces\">\n";
    oss << "     <Topology TopologyType=\"Triangle\" NumberOfElements=\"" << N_Faces << "\">\n";
    oss << "       <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << N_Faces << " 3\">\n";
    oss << "        " << fn.CStr() <<":/FaceCon \n";
    oss << "       </DataItem>\n";
    oss << "     </Topology>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << N_Verts << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Verts \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Cluster\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Cluster \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Float\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngVelocity\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Float\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/AngVelocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    }
    oss << "   <Grid Name=\"DEM_Center\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Particles.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Radius \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Ekin\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PEkin \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PTag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngVel\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PAngVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngMom\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PAngMomentum\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";


    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Domain::WriteFrac (char const * FileKey)
{

    // Counting the number of non valid cohesive interactons
    size_t N_Faces = 0;
    size_t N_Verts = 0;
    size_t nvbi = 0;
    for (size_t i=0;i<BInteractons.Size();i++)
    {
        if (!BInteractons[i]->valid) 
        {
            Particle * P1 = BInteractons[i]->P1;
            Particle * P2 = BInteractons[i]->P2;
            Face     * F1 = P1->Faces[BInteractons[i]->IF1];
            Face     * F2 = P2->Faces[BInteractons[i]->IF2];
            nvbi++;
            N_Faces += F1->Edges.Size();
            N_Faces += F2->Edges.Size();
            N_Verts += F1->Edges.Size() + 1;
            N_Verts += F2->Edges.Size() + 1;
        }
    }

    //std::cout << "1 " << nvbi << std::endl;

    if (nvbi==0) return;

    //Geometric information
    float  * Verts   = new float [3*N_Verts];
    int    * FaceCon = new int   [3*N_Faces];

    //Atributes
    int    * Tags    = new int   [  N_Faces];
       
    size_t n_verts = 0;
    size_t n_faces = 0;
    size_t n_attrs = 0;
    for (size_t i=0;i<BInteractons.Size();i++)
    {
        if (BInteractons[i]->valid) continue;
        Particle * P1 = BInteractons[i]->P1;
        Particle * P2 = BInteractons[i]->P2;
        size_t    IF1 = BInteractons[i]->IF1;
        size_t    IF2 = BInteractons[i]->IF2;
        //std::cout << P1 << " " << P2 <<std::endl;

        //For P1
        {
            size_t n_refv = n_verts/3;
            Array<Vec3_t> Vtemp(P1->Verts.Size());
            Array<Vec3_t> Vres (P1->Verts.Size());
            for (size_t j=0;j<P1->Verts.Size();j++)
            {
                Vtemp[j] = *P1->Verts[j];
                Vres [j] = *P1->Verts[j];
            }
            double multiplier = 0.0;
            if (P1->Eroded&&P1->Faces.Size()>=4)
            {
                DEM::Dilation(Vtemp,P1->EdgeCon,P1->FaceCon,Vres,P1->Props.R);
                multiplier = 1.0;
            }
            for (size_t j=0;j<P1->FaceCon[IF1].Size();j++)
            {
                size_t k = P1->FaceCon[IF1][j];
                Verts[n_verts++] = float(Vres[k](0));
                Verts[n_verts++] = float(Vres[k](1));
                Verts[n_verts++] = float(Vres[k](2));
            }
            size_t n_reff = n_verts/3;
            Vec3_t C,N;
            P1->Faces[IF1]->Centroid(C);
            P1->Faces[IF1]->Normal(N);
            Verts[n_verts++] = float(C(0) + multiplier*P1->Props.R*N(0));
            Verts[n_verts++] = float(C(1) + multiplier*P1->Props.R*N(1));
            Verts[n_verts++] = float(C(2) + multiplier*P1->Props.R*N(2));
            for (size_t j=0;j<P1->FaceCon[IF1].Size();j++)
            {
                FaceCon[n_faces++] = int(n_reff);  
                FaceCon[n_faces++] = int(n_refv + j);
                FaceCon[n_faces++] = int(n_refv + (j+1)%(P1->FaceCon[IF1].Size()));
                Tags   [n_attrs]   = int(P1->Tag);
                n_attrs++;
            }
        }
        //std::cout << "2" << std::endl;
        //For P2
        {
            size_t n_refv = n_verts/3;
            Array<Vec3_t> Vtemp(P2->Verts.Size());
            Array<Vec3_t> Vres (P2->Verts.Size());
            for (size_t j=0;j<P2->Verts.Size();j++)
            {
                Vtemp[j] = *P2->Verts[j];
                Vres [j] = *P2->Verts[j];
            }
            //std::cout << "3" << std::endl;
            double multiplier = 0.0;
            if (P2->Eroded&&P2->Faces.Size()>=4)
            {
                DEM::Dilation(Vtemp,P2->EdgeCon,P2->FaceCon,Vres,P2->Props.R);
                multiplier = 1.0;
            }
            //std::cout << "4" << std::endl;
            for (size_t j=0;j<P2->FaceCon[IF2].Size();j++)
            {
                size_t k = P2->FaceCon[IF2][j];
                Verts[n_verts++] = float(Vres[k](0));
                Verts[n_verts++] = float(Vres[k](1));
                Verts[n_verts++] = float(Vres[k](2));
            }
            //std::cout << "5" << std::endl;
            size_t n_reff = n_verts/3;
            Vec3_t C,N;
            P2->Faces[IF2]->Centroid(C);
            P2->Faces[IF2]->Normal(N);
            Verts[n_verts++] = float(C(0) + multiplier*P2->Props.R*N(0));
            Verts[n_verts++] = float(C(1) + multiplier*P2->Props.R*N(1));
            Verts[n_verts++] = float(C(2) + multiplier*P2->Props.R*N(2));
            //std::cout << "6" << std::endl;
            for (size_t j=0;j<P2->FaceCon[IF2].Size();j++)
            {
                FaceCon[n_faces++] = int(n_reff);  
                FaceCon[n_faces++] = int(n_refv + j);
                FaceCon[n_faces++] = int(n_refv + (j+1)%(P2->FaceCon[IF2].Size()));
                Tags   [n_attrs]   = int(P2->Tag);
                n_attrs++;
            }
            //std::cout << "7" << std::endl;
        }
    }
    //std::cout << n_faces << " " << N_Faces << std::endl;
    //Write the data
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dims[1];
    String dsname;
    dims[0] = 3*N_Verts;
    dsname.Printf("Verts");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Verts);
    dims[0] = 3*N_Faces;
    dsname.Printf("FaceCon");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,FaceCon);
    dims[0] = N_Faces;
    dsname.Printf("Tag");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Tags   );

    //Erasing the data
    delete [] Verts;
    delete [] FaceCon;
    delete [] Tags;

    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"DEM_Faces\">\n";
    oss << "     <Topology TopologyType=\"Triangle\" NumberOfElements=\"" << N_Faces << "\">\n";
    oss << "       <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << N_Faces << " 3\">\n";
    oss << "        " << fn.CStr() <<":/FaceCon \n";
    oss << "       </DataItem>\n";
    oss << "     </Topology>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << N_Verts << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Verts \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";

    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();

}

inline void Domain::Save (char const * FileKey)
{

    // Opening the file for writing
    String fn(FileKey);
    fn.append(".hdf5");
    //if (Util::FileExists(fn))
    //{
        //String command;
        //command.Printf("rm %s",fn.CStr());
        //system(command.CStr());
    //}
    hid_t file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Storing the number of particles in the domain
    int data[1];
    data[0]=Particles.Size();
    hsize_t dims[1];
    dims[0]=1;
    H5LTmake_dataset_int(file_id,"/NP",1,dims,data);

    for (size_t i=0; i<Particles.Size(); i++)
    {
        // Creating the string and the group for each particle
        hid_t group_id;
        String par;
        par.Printf("/Particle_%08d",i);
        group_id = H5Gcreate(file_id, par.CStr(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


        // Storing some scalar variables
        double dat[1];
        dat[0] = Particles[i]->Props.R;
        H5LTmake_dataset_double(group_id,"SR",1,dims,dat);
        dat[0] = Particles[i]->Props.rho;
        H5LTmake_dataset_double(group_id,"Rho",1,dims,dat);
        dat[0] = Particles[i]->Props.m;
        H5LTmake_dataset_double(group_id,"m",1,dims,dat);
        dat[0] = Particles[i]->Props.V;
        H5LTmake_dataset_double(group_id,"V",1,dims,dat);
        dat[0] = Particles[i]->Diam;
        H5LTmake_dataset_double(group_id,"Diam",1,dims,dat);
        dat[0] = Particles[i]->Dmax;
        H5LTmake_dataset_double(group_id,"Dmax",1,dims,dat);
        int datint[1];
        datint[0] = Particles[i]->Index;
        H5LTmake_dataset_int(group_id,"Index",1,dims,datint);


        int tag[1];
        tag[0] = Particles[i]->Tag;
        H5LTmake_dataset_int(group_id,"Tag",1,dims,tag);

        // Storing vectorial variables
        double cd[3];
        hsize_t dd[1];
        dd[0] = 3;

        cd[0]=Particles[i]->x(0);
        cd[1]=Particles[i]->x(1);
        cd[2]=Particles[i]->x(2);
        H5LTmake_dataset_double(group_id,"x",1,dd,cd);

        cd[0]=Particles[i]->xb(0);
        cd[1]=Particles[i]->xb(1);
        cd[2]=Particles[i]->xb(2);
        H5LTmake_dataset_double(group_id,"xb",1,dd,cd);

        cd[0]=Particles[i]->v(0);
        cd[1]=Particles[i]->v(1);
        cd[2]=Particles[i]->v(2);
        H5LTmake_dataset_double(group_id,"v",1,dd,cd);

        cd[0]=Particles[i]->w(0);
        cd[1]=Particles[i]->w(1);
        cd[2]=Particles[i]->w(2);
        H5LTmake_dataset_double(group_id,"w",1,dd,cd);

        cd[0]=Particles[i]->wb(0);
        cd[1]=Particles[i]->wb(1);
        cd[2]=Particles[i]->wb(2);
        H5LTmake_dataset_double(group_id,"wb",1,dd,cd);

        cd[0]=Particles[i]->I(0);
        cd[1]=Particles[i]->I(1);
        cd[2]=Particles[i]->I(2);
        H5LTmake_dataset_double(group_id,"I",1,dd,cd);

        double cq[4];
        dd[0] = 4;
        cq[0]=Particles[i]->Q(0);
        cq[1]=Particles[i]->Q(1);
        cq[2]=Particles[i]->Q(2);
        cq[3]=Particles[i]->Q(3);
        H5LTmake_dataset_double(group_id,"Q",1,dd,cq);




        // Storing the number of vertices of each particle
        data[0] = Particles[i]->Verts.Size();
        H5LTmake_dataset_int(group_id,"n_vertices",1,dims,data);
        hid_t gv_id;
        gv_id = H5Gcreate(group_id,"Verts", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Storing each vertex 
        for (size_t j=0;j<Particles[i]->Verts.Size();j++)
        {
            String parv;
            parv.Printf("Verts_%08d",j);
            double cod[3];
            cod[0]=(*Particles[i]->Verts[j])(0);
            cod[1]=(*Particles[i]->Verts[j])(1);
            cod[2]=(*Particles[i]->Verts[j])(2);
            hsize_t dim[1];
            dim[0]=3;
            H5LTmake_dataset_double(gv_id,parv.CStr(),1,dim,cod);
        }

        // Number of edges of the particle
        data[0] = Particles[i]->Edges.Size();
        H5LTmake_dataset_int(group_id,"n_edges",1,dims,data);
        gv_id = H5Gcreate(group_id,"Edges", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Edges
        for (size_t j=0;j<Particles[i]->Edges.Size();j++)
        {
            String parv;
            parv.Printf("Edges_%08d",j);
            int co[2];
            co[0] = Particles[i]->EdgeCon[j][0];
            co[1] = Particles[i]->EdgeCon[j][1];
            hsize_t dim[1];
            dim[0] =2;
            H5LTmake_dataset_int(gv_id,parv.CStr(),1,dim,co);
        }
        
        // Number of faces of the particle
        data[0] = Particles[i]->Faces.Size();
        H5LTmake_dataset_int(group_id,"n_faces",1,dims,data);
        gv_id = H5Gcreate(group_id,"Faces", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Number of cylinders of the particle
        data[0] = Particles[i]->Cylinders.Size();
        H5LTmake_dataset_int(group_id,"n_cylinders",1,dims,data);
        
        // Faces
        for (size_t j=0;j<Particles[i]->Faces.Size();j++)
        {
            String parv;
            parv.Printf("Faces_%08d",j);
            int co[Particles[i]->FaceCon[j].Size()];
            hsize_t dim[1];
            dim[0]= Particles[i]->FaceCon[j].Size();
            for (size_t k=0;k<Particles[i]->FaceCon[j].Size();k++)
            {
                co[k]=Particles[i]->FaceCon[j][k];
            }
            H5LTmake_dataset_int(gv_id,parv.CStr(),1,dim,co);
        }
        
    }

    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);
    //sleep(5);
}

inline void Domain::Load (char const * FileKey)
{

    // Opening the file for reading
    String fn(FileKey);
    fn.append(".hdf5");
    if (!Util::FileExists(fn)) throw new Fatal("File <%s> not found",fn.CStr());
    printf("\n%s--- Loading file %s --------------------------------------------%s\n",TERM_CLR1,fn.CStr(),TERM_RST);
    hid_t file_id;
    file_id = H5Fopen(fn.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // Number of particles in the domain
    int data[1];
    H5LTread_dataset_int(file_id,"/NP",data);
    size_t NP = data[0];

    // Loading the particles
    for (size_t i=0; i<NP; i++)
    {

        // Creating the string and the group for each particle
        hid_t group_id;
        String par;
        par.Printf("/Particle_%08d",i);
        group_id = H5Gopen(file_id, par.CStr(),H5P_DEFAULT);

        // Finding the particle's position for the domain decomposition
        double X[3];
        H5LTread_dataset_double(group_id,"x",X);


        // Loading the Vertices
        H5LTread_dataset_int(group_id,"n_vertices",data);
        size_t nv = data[0];
        hid_t gv_id;
        gv_id = H5Gopen(group_id,"Verts", H5P_DEFAULT);
        Array<Vec3_t> V;

        for (size_t j=0;j<nv;j++)
        {
            String parv;
            parv.Printf("Verts_%08d",j);
            double cod[3];
            H5LTread_dataset_double(gv_id,parv.CStr(),cod);
            V.Push(Vec3_t(cod[0],cod[1],cod[2]));
        }
        
        // Loading the edges
        H5LTread_dataset_int(group_id,"n_edges",data);
        size_t ne = data[0];
        gv_id = H5Gopen(group_id,"Edges", H5P_DEFAULT);
        Array<Array <int> > E;

        for (size_t j=0;j<ne;j++)
        {
            String parv;
            parv.Printf("Edges_%08d",j);
            int cod[2];
            H5LTread_dataset_int(gv_id,parv.CStr(),cod);
            Array<int> Ep(2);
            Ep[0]=cod[0];
            Ep[1]=cod[1];
            E.Push(Ep);
        }

        // Loading the faces

        // Number of faces of the particle
        H5LTread_dataset_int(group_id,"n_faces",data);
        size_t nf = data[0];
        gv_id = H5Gopen(group_id,"Faces", H5P_DEFAULT);
        Array<Array <int> > F;
        
        // Faces
        for (size_t j=0;j<nf;j++)
        {
            String parv;
            parv.Printf("Faces_%08d",j);
            hsize_t dim[1];
            H5LTget_dataset_info(gv_id,parv.CStr(),dim,NULL,NULL);
            size_t ns = (size_t)dim[0];
            int co[ns];
            Array<int> Fp(ns);

            H5LTread_dataset_int(gv_id,parv.CStr(),co);
            
            for (size_t k=0;k<ns;k++)
            {
                Fp[k] = co[k];
            }

            F.Push(Fp);

        }

        // Number of cylinders
        H5LTread_dataset_int(group_id,"n_cylinders",data);
        size_t nc = data[0];

        Particles.Push (new Particle(-1,V,E,F,OrthoSys::O,OrthoSys::O,0.1,1.0));

        // Loading cylinder data if applicable
        if (nc>0)
        {   
            Vec3_t X0 = 0.5*(*Particles[Particles.Size()-1]->Verts[0] + *Particles[Particles.Size()-1]->Verts[2]);
            Vec3_t X1 = 0.5*(*Particles[Particles.Size()-1]->Verts[3] + *Particles[Particles.Size()-1]->Verts[5]);
            Particles[Particles.Size()-1]->Tori.Push     (new Torus(&X0,Particles[Particles.Size()-1]->Verts[0],Particles[Particles.Size()-1]->Verts[1]));
            Particles[Particles.Size()-1]->Tori.Push     (new Torus(&X1,Particles[Particles.Size()-1]->Verts[3],Particles[Particles.Size()-1]->Verts[4]));
            Particles[Particles.Size()-1]->Cylinders.Push(new Cylinder(Particles[Particles.Size()-1]->Tori[0],Particles[Particles.Size()-1]->Tori[1],Particles[Particles.Size()-1]->Verts[2],Particles[Particles.Size()-1]->Verts[5]));
        }

        // Loading vectorial variables
        Particles[Particles.Size()-1]->x = Vec3_t(X[0],X[1],X[2]);
        double cd[3];
        H5LTread_dataset_double(group_id,"xb",cd);
        Particles[Particles.Size()-1]->xb = Vec3_t(cd[0],cd[1],cd[2]);
        H5LTread_dataset_double(group_id,"v",cd);
        Particles[Particles.Size()-1]->v = Vec3_t(cd[0],cd[1],cd[2]);
        H5LTread_dataset_double(group_id,"w",cd);
        Particles[Particles.Size()-1]->w = Vec3_t(cd[0],cd[1],cd[2]);
        H5LTread_dataset_double(group_id,"wb",cd);
        Particles[Particles.Size()-1]->wb = Vec3_t(cd[0],cd[1],cd[2]);
        H5LTread_dataset_double(group_id,"I",cd);
        Particles[Particles.Size()-1]->I = Vec3_t(cd[0],cd[1],cd[2]);

        double cq[4];
        H5LTread_dataset_double(group_id,"Q",cq);
        Particles[Particles.Size()-1]->Q = Quaternion_t(cq[0],cq[1],cq[2],cq[3]);
    
        // Loading the scalar quantities of the particle
        double dat[1];
        H5LTread_dataset_double(group_id,"SR",dat);
        Particles[Particles.Size()-1]->Props.R = dat[0];
        H5LTread_dataset_double(group_id,"Rho",dat);
        Particles[Particles.Size()-1]->Props.rho = dat[0];
        H5LTread_dataset_double(group_id,"m",dat);
        Particles[Particles.Size()-1]->Props.m = dat[0];
        H5LTread_dataset_double(group_id,"V",dat);
        Particles[Particles.Size()-1]->Props.V = dat[0];
        H5LTread_dataset_double(group_id,"Diam",dat);
        Particles[Particles.Size()-1]->Diam = dat[0];
        H5LTread_dataset_double(group_id,"Dmax",dat);
        Particles[Particles.Size()-1]->Dmax = dat[0];
        int datint[1];
        H5LTread_dataset_int(group_id,"Index",datint);
        //Particles[Particles.Size()-1]->Index = datint[0];
        Particles[Particles.Size()-1]->Index = Particles.Size()-1;
        int tag[1];
        H5LTread_dataset_int(group_id,"Tag",tag);
        Particles[Particles.Size()-1]->Tag = tag[0];
        Particles[Particles.Size()-1]->PropsReady = true;

    }


    H5Fclose(file_id);
    printf("\n%s--- Done --------------------------------------------%s\n",TERM_CLR2,TERM_RST);
}

#endif

#ifdef USE_VTK

void Domain::WriteVTKContacts  (char const * FileKey)
{
    size_t ncontacts = 0;

    for (size_t i=0;i<CInteractons.Size();i++)
    {
        if ((norm(CInteractons[i]->Fnet)>0.0)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree())) ncontacts++;
    }

    if (ncontacts==0) return;
    
    // Create a vtkPoints object and store the points in it
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToDouble(); // set the precision to double
    // add the points     // WARNING: not considering periodic conditions or walls
    for (size_t i=0;i<CInteractons.Size();i++)
    {
        if ((norm(CInteractons[i]->Fnet)>0.0)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree()))
        {
            Vec3_t R1 = CInteractons[i]->P1->x; 
            Vec3_t R2 = CInteractons[i]->P2->x; 
            points->InsertNextPoint(R1(0), R1(1), R1(2)); 
            points->InsertNextPoint(R2(0), R2(1), R2(2)); 
        }
    } 
    const auto npoints = points->GetNumberOfPoints();
    
    // Create a cell array to store the lines in and add the lines to it
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    // create each line according to contact
    for(int ip = 0; ip < npoints; ip += 2) 
    {
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, ip);
        line->GetPointIds()->SetId(1, ip+1);
        lines->InsertNextCell(line);
    }
    
    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    // Add the points to the dataset
    polyData->SetPoints(points);
    // Add the lines to the dataset
    polyData->SetLines(lines);
    
    // ADD CELL DATA
    // SEE: http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/MiscCellData
    // SCALARS
    /*// uid1
    vtkSmartPointer<vtkIntArray> uid1 = 
      vtkSmartPointer<vtkIntArray>::New();
    uid1->SetNumberOfComponents(1); 
    uid1->SetNumberOfTuples(ncontacts); 
    uid1->SetName("uid1");
    for (int ic = 0; ic < ncontacts; ++ic) {
      uid1->SetValue(ic, contacts[ic].uid1_); 
    }
    polyData->GetCellData()->AddArray(uid1);*/
    // VECTORS
    // Normal
    vtkSmartPointer<vtkDoubleArray> NormalForce = 
    vtkSmartPointer<vtkDoubleArray>::New();
    NormalForce->SetNumberOfComponents(1); 
    NormalForce->SetNumberOfTuples(ncontacts); 
    NormalForce->SetName("Normal Force");
    size_t ic = 0;
    for (size_t i=0;i<CInteractons.Size();i++)
    {
        if ((norm(CInteractons[i]->Fnet)>0.0)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree()))
        {
            Vec3_t P = CInteractons[i]->Fnet;
            double data[1] = {norm(P)}; 
            NormalForce->SetTupleValue(ic, data); 
            ic++;
        }
    }
    polyData->GetCellData()->AddArray(NormalForce);
    
    // Write the file
    String fn(FileKey);
    fn.append(".vtp");
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(fn.CStr());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(polyData);
#else
    writer->SetInputData(polydata);
#endif
    //// Optional - set the mode.
    writer->SetDataModeToBinary(); // default
    //writer->SetDataModeToAscii();
    writer->SetCompressorTypeToZLib(); // or ToNone
    // write
    writer->Write();    
}

#endif // USE_VTK

#ifdef USE_OMP

inline void Domain::UpdateLinkedCells()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].LPP.Resize(0);
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<FreePar.Size()  ;i++)
    for (size_t j=0;j<NoFreePar.Size();j++)
    {
        size_t i1 = std::min(FreePar[i],NoFreePar[j]);
        size_t i2 = std::max(FreePar[i],NoFreePar[j]);
        MTD[omp_get_thread_num()].LPP.Push(std::make_pair(i1,i2));
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t idx=0;idx<LinkedCell.Size();idx++)
    {
        if (LinkedCell[idx].Size()==0) continue;
        iVec3_t index;
        idx2Pt(idx,index,LCellDim);
        //std::cout << index << " " << LinkedCell[idx].Size() << " ";
        for (size_t n=0  ;n<LinkedCell[idx].Size()-1;n++)
        {
            //std::cout << LinkedCell[idx][n] << " ";
            for (size_t m=n+1;m<LinkedCell[idx].Size()  ;m++)
            {
                size_t i1 = LinkedCell[idx][n];
                size_t i2 = LinkedCell[idx][m];
                if (i1==i2) continue;
                MTD[omp_get_thread_num()].LPP.Push(std::make_pair(i1,i2));
            }
        }
        //std::cout << std::endl;
        size_t i = index(0);
        size_t j = index(1);
        size_t k = index(2);
        for (size_t knb=std::max(0,int(k)-1);knb<=std::min(LCellDim(2)-1,k+1);knb++)
        for (size_t jnb=std::max(0,int(j)-1);jnb<=std::min(LCellDim(1)-1,j+1);jnb++)
        for (size_t inb=std::max(0,int(i)-1);inb<=std::min(LCellDim(0)-1,i+1);inb++)
        {
            iVec3_t Ptnb(inb,jnb,knb);
            size_t idxnb = Pt2idx(Ptnb,LCellDim);
            if (idxnb>idx)
            {
                for (size_t n=0;n<LinkedCell[idx].Size()  ;n++)
                {
                    for (size_t m=0;m<LinkedCell[idxnb].Size()  ;m++)
                    {
                        size_t i1 = std::min(LinkedCell[idx  ][n],LinkedCell[idxnb][m]);
                        size_t i2 = std::max(LinkedCell[idx  ][n],LinkedCell[idxnb][m]);
                        if (i1==i2) continue;
                        MTD[omp_get_thread_num()].LPP.Push(std::make_pair(i1,i2));
                    }
                }
            }
        }
    }
    size_t Npp = 0;
    for (size_t i=0;i<Nproc;i++)
    {
        Npp += MTD[i].LPP.Size();
    }
    ListPosPairs.Resize(Npp);
    size_t idx = 0;
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LPP.Size();j++)
        {
            ListPosPairs[idx] = MTD[i].LPP[j];
            idx++;
        }
    }
}

#endif

inline void Domain::BoundingBox(Vec3_t & minX, Vec3_t & maxX)
{
    if (Particles.Size()==0) throw new Fatal("DEM::Domain::BoundingBox: There are no particles to build the bounding box");
    minX = Vec3_t(Particles[0]->MinX(), Particles[0]->MinY(), Particles[0]->MinZ());
    maxX = Vec3_t(Particles[0]->MaxX(), Particles[0]->MaxY(), Particles[0]->MaxZ());
    for (size_t i=1; i<Particles.Size(); i++)
    {
        if (minX(0)>Particles[i]->MinX()&&Particles[i]->IsFree()) minX(0) = Particles[i]->MinX();
        if (minX(1)>Particles[i]->MinY()&&Particles[i]->IsFree()) minX(1) = Particles[i]->MinY();
        if (minX(2)>Particles[i]->MinZ()&&Particles[i]->IsFree()) minX(2) = Particles[i]->MinZ();
        if (maxX(0)<Particles[i]->MaxX()&&Particles[i]->IsFree()) maxX(0) = Particles[i]->MaxX();
        if (maxX(1)<Particles[i]->MaxY()&&Particles[i]->IsFree()) maxX(1) = Particles[i]->MaxY();
        if (maxX(2)<Particles[i]->MaxZ()&&Particles[i]->IsFree()) maxX(2) = Particles[i]->MaxZ();
    }
}

inline void Domain::Center(Vec3_t C)
{
    Vec3_t minX,maxX;
    BoundingBox(minX,maxX);
    Vec3_t Transport(-0.5*(maxX+minX));
    Transport += C;
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Translate(Transport);
}

inline void Domain::ClearInteractons()
{
    // delete old interactors
    for (size_t i=0; i<CInteractons.Size(); ++i)
    {
        if (CInteractons[i]!=NULL) delete CInteractons[i];
    }
    CInteractons.Resize(0);
    for (size_t i=0; i<BInteractons.Size(); ++i)
    {
        if (BInteractons[i]!=NULL) delete BInteractons[i];
    }
    BInteractons.Resize(0);
    Interactons.Resize(0);
    Listofpairs.clear();

    for (size_t i=0; i<CPxInteractons.Size(); ++i)
    {
        if (CPxInteractons[i]!=NULL) delete CPxInteractons[i];
    }
    CPxInteractons.Resize(0);
    PxInteractons.Resize(0);
    PxListofpairs.clear();
    for (size_t i=0; i<CPyInteractons.Size(); ++i)
    {
        if (CPyInteractons[i]!=NULL) delete CPyInteractons[i];
    }
    CPyInteractons.Resize(0);
    PyInteractons.Resize(0);
    PyListofpairs.clear();
}

inline void Domain::ResetInteractons()
{
    // delete old interactors
    for (size_t i=0; i<CInteractons.Size(); ++i)
    {
        if (CInteractons[i]!=NULL) delete CInteractons[i];
    }

    // new interactors
    CInteractons.Resize(0);
    for (size_t i=0; i<Particles.Size()-1; i++)
    {
        bool pi_has_vf = !Particles[i]->IsFree();
        for (size_t j=i+1; j<Particles.Size(); j++)
        {
            bool pj_has_vf = !Particles[j]->IsFree();


            // if both particles have any component specified or they are far away, don't create any intereactor
            bool close = (Distance(Particles[i]->x,Particles[j]->x)<=Particles[i]->Dmax+Particles[j]->Dmax+2*Alpha);
            if ((pi_has_vf && pj_has_vf) || !close ) continue;
            Listofpairs.insert(std::make_pair(Particles[i],Particles[j]));

            // if both particles are spheres (just one vertex)
            if (Particles[i]->Verts.Size()==1 && Particles[j]->Verts.Size()==1)
            {
                CInteractons.Push (new CInteractonSphere(Particles[i],Particles[j]));
            }

            // normal particles
            else
            {
                CInteractons.Push (new CInteracton(Particles[i],Particles[j]));
            }
        }
    }
}

inline void Domain::ResetDisplacements()
{
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].Dmx = 0.0;
        MTD[i].LLC.Resize(0);
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Particles.Size();i++)
    {
        Particles[i]->ResetDisplacements();
        if(Particles[i]->IsFree())
        {
            iVec3_t idx = (Particles[i]->x - LCxmin)/(2.0*Beta*MaxDmax);
            MTD[omp_get_thread_num()].LLC.Push(std::make_pair(idx,i));
        }
    }
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LLC.Size();j++)
        {
            size_t idx = Pt2idx(MTD[i].LLC[j].first,LCellDim);
            LinkedCell[idx].Push(MTD[i].LLC[j].second);
        }
    }
#else
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Particles[i]->ResetDisplacements();
    }
#endif
}

inline double Domain::MaxDisplacement()
{
    double md = 0.0;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        double mpd = Particles[i]->MaxDisplacement();
        if (mpd > md) md = mpd;
    }
    return md;
}

inline void Domain::ResetContacts()
{   
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].LC.Resize(0);
        MTD[i].LCI.Resize(0);
        MTD[i].LCB.Resize(0);
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t n=0;n<ListPosPairs.Size();n++)
    {
        size_t i = ListPosPairs[n].first;
        size_t j = ListPosPairs[n].second;
        bool pi_has_vf = !Particles[i]->IsFree();
        bool pj_has_vf = !Particles[j]->IsFree();
        bool close = (Distance(Particles[i]->x,Particles[j]->x)<=Particles[i]->Dmax+Particles[j]->Dmax+2*Alpha);
        if ((pi_has_vf && pj_has_vf) || !close) continue;
        std::set<std::pair<Particle *, Particle *> >::iterator it = Listofpairs.find(std::make_pair(Particles[i],Particles[j]));
        if (it != Listofpairs.end())
        {
            continue;
        }
        MTD[omp_get_thread_num()].LC.Push(std::make_pair(i,j));
    }
    for (size_t i=0;i<Nproc;i++)
    {
        //std::cout << MTD[i].LC.Size() << std::endl;
        for (size_t j=0;j<MTD[i].LC.Size();j++)
        {
        //std::cout << MTD[i].LC.Size() << std::endl;
            size_t n = MTD[i].LC[j].first;
            size_t m = MTD[i].LC[j].second;
            Listofpairs.insert(std::make_pair(Particles[n],Particles[m]));
            if (Particles[n]->Verts.Size()==1 && Particles[m]->Verts.Size()==1)
            {
                if (!MostlySpheres) CInteractons.Push (new CInteractonSphere(Particles[n],Particles[m]));
                //if (!MostlySpheres) SInteractons.Push (new SphereCollision(Particles[n],Particles[m]));
            }
            else
            {
                CInteractons.Push (new CInteracton(Particles[n],Particles[m]));
            }
        }
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t n=0;n<CInteractons.Size();n++)
    {
        if(CInteractons[n]->UpdateContacts(Alpha)) MTD[omp_get_thread_num()].LCI.Push(n);
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t n=0;n<BInteractons.Size();n++)
    {
        if(BInteractons[n]->UpdateContacts(Alpha)) MTD[omp_get_thread_num()].LCB.Push(n);
    }
    Interactons.Resize(0);
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LCI.Size();j++)
        {
            Interactons.Push(CInteractons[MTD[i].LCI[j]]);
        }
        for (size_t j=0;j<MTD[i].LCB.Size();j++)
        {
            Interactons.Push(BInteractons[MTD[i].LCB[j]]);
        }
    }
    if (Xmax-Xmin>Alpha)
    {
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<Nproc;i++)
        {
            MTD[i].LBP.Resize(0);
            MTD[i].LPC.Resize(0);
            MTD[i].LPCI.Resize(0);
        }

        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<Particles.Size();i++)
        {
            Particle * Pa = Particles[i];
            if ((Pa->MinX()>Xmax)&&Pa->IsFree())
            {
                Vec3_t v(Xmin-Xmax,0.0,0.0);
                Pa->Translate(v);
            }
            if ((Pa->MaxX()<Xmin)&&Pa->IsFree())
            {
                Vec3_t v(Xmax-Xmin,0.0,0.0);
                Pa->Translate(v);
            }
            if ((Pa->MaxX()>Xmax-2.0*Alpha-2.0*MaxDmax)&&Pa->IsFree())
            {
                MTD[omp_get_thread_num()].LBP.Push(i);
            }
        }
        ParXmax.Resize(0);
        for (size_t i=0;i<Nproc;i++)
        {
            for (size_t j=0;j<MTD[i].LBP.Size();j++)
            {
                ParXmax.Push(Particles[MTD[i].LBP[j]]);
            }
        }

        Vec3_t v(Xmin-Xmax,0.0,0.0);
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<ParXmax.Size();i++)
        {
            Particle * P1 = ParXmax[i];
            P1->Translate(v);
            iVec3_t index = (P1->x - LCxmin)/(2.0*Beta*MaxDmax);
            size_t ic = index(0);
            size_t jc = index(1);
            size_t kc = index(2);
            if (index(0)>LCellDim(0)) ic = 0;
            for (size_t knb=std::max(0,int(kc)-1);knb<=std::min(LCellDim(2)-1,kc+1);knb++)
            for (size_t jnb=std::max(0,int(jc)-1);jnb<=std::min(LCellDim(1)-1,jc+1);jnb++)
            for (size_t inb=std::max(0,int(ic)-1);inb<=std::min(LCellDim(0)-1,ic+1);inb++)
            {
                iVec3_t Ptnb(inb,jnb,knb);
                size_t idxnb = Pt2idx(Ptnb,LCellDim);
                for (size_t m=0;m<LinkedCell[idxnb].Size()  ;m++)
                {
                    size_t i2 = LinkedCell[idxnb][m];
                    Particle * P2 = Particles[i2];
                    if (P1==P2||ParXmax.Has(P2))
                    {
                        //std::cout << P1->x << " " << P2->x << std::endl;
                        continue;
                    }

                    bool close = (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*Alpha);
                    if (!close) continue;
                    std::set<std::pair<Particle *, Particle *> >::iterator it = PxListofpairs.find(std::make_pair(P1,P2));
                    if (it != PxListofpairs.end())
                    {
                        continue;
                    }
                    MTD[omp_get_thread_num()].LPC.Push(std::make_pair(i,i2));
                }
            }
        }

        for (size_t i=0;i<Nproc;i++)
        {
            for (size_t j=0;j<MTD[i].LPC.Size();j++)
            {
                size_t n = MTD[i].LPC[j].first;
                size_t m = MTD[i].LPC[j].second;
                PxListofpairs.insert(std::make_pair(ParXmax[n],Particles[m]));
                if (ParXmax[n]->Verts.Size()==1 && Particles[m]->Verts.Size()==1)
                {
                    CPxInteractons.Push (new CInteractonSphere(ParXmax[n],Particles[m]));
                }
                else
                {
                    CPxInteractons.Push (new CInteracton(ParXmax[n],Particles[m]));
                }
            }
        }

        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t n=0;n<CPxInteractons.Size();n++)
        {
            if(CPxInteractons[n]->UpdateContacts(Alpha)) MTD[omp_get_thread_num()].LPCI.Push(n);
        }
        PxInteractons.Resize(0);
        for (size_t i=0;i<Nproc;i++)
        {
            for (size_t j=0;j<MTD[i].LPCI.Size();j++)
            {
                PxInteractons.Push(CPxInteractons[MTD[i].LPCI[j]]);
            }
        }
        v = Vec3_t(Xmax-Xmin,0.0,0.0);
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<ParXmax.Size();i++)
        {
            Particle * P1 = ParXmax[i];
            P1->Translate(v);
        }
    }
    if (Ymax-Ymin>Alpha)
    {
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<Nproc;i++)
        {
            MTD[i].LBP.Resize(0);
            MTD[i].LPC.Resize(0);
            MTD[i].LPCI.Resize(0);
        }

        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<Particles.Size();i++)
        {
            Particle * Pa = Particles[i];
            if ((Pa->MinY()>Ymax)&&Pa->IsFree())
            {
                Vec3_t v(0.0,Ymin-Ymax,0.0);
                Pa->Translate(v);
            }
            if ((Pa->MaxY()<Ymin)&&Pa->IsFree())
            {
                Vec3_t v(0.0,Ymax-Ymin,0.0);
                Pa->Translate(v);
            }
            if ((Pa->MaxY()>Ymax-2.0*Alpha-2.0*MaxDmax)&&Pa->IsFree())
            {
                MTD[omp_get_thread_num()].LBP.Push(i);
            }
        }
        ParYmax.Resize(0);
        for (size_t i=0;i<Nproc;i++)
        {
            for (size_t j=0;j<MTD[i].LBP.Size();j++)
            {
                ParYmax.Push(Particles[MTD[i].LBP[j]]);
            }
        }

        Vec3_t v(0.0,Ymin-Ymax,0.0);
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<ParYmax.Size();i++)
        {
            Particle * P1 = ParYmax[i];
            P1->Translate(v);
            iVec3_t index = (P1->x - LCxmin)/(2.0*Beta*MaxDmax);
            size_t ic = index(0);
            size_t jc = index(1);
            size_t kc = index(2);
            if (index(1)>LCellDim(1)) jc = 0;
            for (size_t knb=std::max(0,int(kc)-1);knb<=std::min(LCellDim(2)-1,kc+1);knb++)
            for (size_t jnb=std::max(0,int(jc)-1);jnb<=std::min(LCellDim(1)-1,jc+1);jnb++)
            for (size_t inb=std::max(0,int(ic)-1);inb<=std::min(LCellDim(0)-1,ic+1);inb++)
            {
                iVec3_t Ptnb(inb,jnb,knb);
                size_t idxnb = Pt2idx(Ptnb,LCellDim);
                for (size_t m=0;m<LinkedCell[idxnb].Size()  ;m++)
                {
                    size_t i2 = LinkedCell[idxnb][m];
                    Particle * P2 = Particles[i2];
                    if (P1==P2||ParYmax.Has(P2))
                    {
                        //std::cout << P1->x << " " << P2->x << std::endl;
                        continue;
                    }

                    bool close = (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*Alpha);
                    if (!close) continue;
                    std::set<std::pair<Particle *, Particle *> >::iterator it = PyListofpairs.find(std::make_pair(P1,P2));
                    if (it != PyListofpairs.end())
                    {
                        continue;
                    }
                    MTD[omp_get_thread_num()].LPC.Push(std::make_pair(i,i2));
                }
            }
        }

        for (size_t i=0;i<Nproc;i++)
        {
            for (size_t j=0;j<MTD[i].LPC.Size();j++)
            {
                size_t n = MTD[i].LPC[j].first;
                size_t m = MTD[i].LPC[j].second;
                PyListofpairs.insert(std::make_pair(ParYmax[n],Particles[m]));
                if (ParYmax[n]->Verts.Size()==1 && Particles[m]->Verts.Size()==1)
                {
                    CPyInteractons.Push (new CInteractonSphere(ParYmax[n],Particles[m]));
                }
                else
                {
                    CPyInteractons.Push (new CInteracton(ParYmax[n],Particles[m]));
                }
            }
        }

        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t n=0;n<CPyInteractons.Size();n++)
        {
            if(CPyInteractons[n]->UpdateContacts(Alpha)) MTD[omp_get_thread_num()].LPCI.Push(n);
        }
        PyInteractons.Resize(0);
        for (size_t i=0;i<Nproc;i++)
        {
            for (size_t j=0;j<MTD[i].LPCI.Size();j++)
            {
                PyInteractons.Push(CPyInteractons[MTD[i].LPCI[j]]);
            }
        }
        v = Vec3_t(0.0,Ymax-Ymin,0.0);
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<ParYmax.Size();i++)
        {
            Particle * P1 = ParYmax[i];
            P1->Translate(v);
        }
    }
    if ((Xmax-Xmin>Alpha)&&(Ymax-Ymin>Alpha))
    {
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<Nproc;i++)
        {
            MTD[i].LBP.Resize(0);
            MTD[i].LPC.Resize(0);
            MTD[i].LPCI.Resize(0);
        }

        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<ParXmax.Size();i++)
        {
            if (ParYmax.Has(ParXmax[i]))
            {
                MTD[omp_get_thread_num()].LBP.Push(ParXmax[i]->Index);
            }
        }

        ParXYmax.Resize(0);
        for (size_t i=0;i<Nproc;i++)
        {
            for (size_t j=0;j<MTD[i].LBP.Size();j++)
            {
                ParXYmax.Push(Particles[MTD[i].LBP[j]]);
            }
        }
        Vec3_t v(Xmin-Xmax,Ymin-Ymax,0.0);
        for (size_t i=0;i<ParXYmax.Size();i++)
        {
            Particle * P1 = ParXYmax[i];
            P1->Translate(v);
            iVec3_t index = (P1->x - LCxmin)/(2.0*Beta*MaxDmax);
            size_t ic = index(0);
            size_t jc = index(1);
            size_t kc = index(2);
            if (index(0)>LCellDim(0)) ic = 0;
            if (index(1)>LCellDim(1)) jc = 0;
            for (size_t knb=std::max(0,int(kc)-1);knb<=std::min(LCellDim(2)-1,kc+1);knb++)
            for (size_t jnb=std::max(0,int(jc)-1);jnb<=std::min(LCellDim(1)-1,jc+1);jnb++)
            for (size_t inb=std::max(0,int(ic)-1);inb<=std::min(LCellDim(0)-1,ic+1);inb++)
            {
                iVec3_t Ptnb(inb,jnb,knb);
                size_t idxnb = Pt2idx(Ptnb,LCellDim);
                for (size_t m=0;m<LinkedCell[idxnb].Size()  ;m++)
                {
                    size_t i2 = LinkedCell[idxnb][m];
                    Particle * P2 = Particles[i2];
                    if (P1==P2||ParXYmax.Has(P2))
                    {
                        //std::cout << P1->x << " " << P2->x << std::endl;
                        continue;
                    }

                    bool close = (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*Alpha);
                    if (!close) continue;
                    std::set<std::pair<Particle *, Particle *> >::iterator it = PxyListofpairs.find(std::make_pair(P1,P2));
                    if (it != PxyListofpairs.end())
                    {
                        continue;
                    }
                    MTD[omp_get_thread_num()].LPC.Push(std::make_pair(i,i2));
                }
            }
        }
        for (size_t i=0;i<Nproc;i++)
        {
            for (size_t j=0;j<MTD[i].LPC.Size();j++)
            {
                size_t n = MTD[i].LPC[j].first;
                size_t m = MTD[i].LPC[j].second;
                PxyListofpairs.insert(std::make_pair(ParXYmax[n],Particles[m]));
                if (ParXYmax[n]->Verts.Size()==1 && Particles[m]->Verts.Size()==1)
                {
                    CPxyInteractons.Push (new CInteractonSphere(ParXYmax[n],Particles[m]));
                }
                else
                {
                    CPxyInteractons.Push (new CInteracton(ParXYmax[n],Particles[m]));
                }
            }
        }

        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t n=0;n<CPxyInteractons.Size();n++)
        {
            if(CPxyInteractons[n]->UpdateContacts(Alpha)) MTD[omp_get_thread_num()].LPCI.Push(n);
        }
        PxyInteractons.Resize(0);
        for (size_t i=0;i<Nproc;i++)
        {
            for (size_t j=0;j<MTD[i].LPCI.Size();j++)
            {
                PxyInteractons.Push(CPxyInteractons[MTD[i].LPCI[j]]);
            }
        }
        v = Vec3_t(Xmax-Xmin,Ymax-Ymin,0.0);
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<ParXYmax.Size();i++)
        {
            Particle * P1 = ParXYmax[i];
            P1->Translate(v);
        }
    }
#else
    for (size_t i=0; i<Particles.Size()-1; i++)
    {
        bool pi_has_vf = !Particles[i]->IsFree();
        for (size_t j=i+1; j<Particles.Size(); j++)
        {
            bool pj_has_vf = !Particles[j]->IsFree();

            bool close = (Distance(Particles[i]->x,Particles[j]->x)<=Particles[i]->Dmax+Particles[j]->Dmax+2*Alpha);
            if ((pi_has_vf && pj_has_vf) || !close) continue;
            
            // checking if the interacton exist for that pair of particles
            std::set<std::pair<Particle *, Particle *> >::iterator it = Listofpairs.find(std::make_pair(Particles[i],Particles[j]));
            if (it != Listofpairs.end())
            {
                continue;
            }
            Listofpairs.insert(std::make_pair(Particles[i],Particles[j]));
            
            // if both particles are spheres (just one vertex)
            if (Particles[i]->Verts.Size()==1 && Particles[j]->Verts.Size()==1)
            {
                CInteractons.Push (new CInteractonSphere(Particles[i],Particles[j]));
            }

            // normal particles
            else
            {
                CInteractons.Push (new CInteracton(Particles[i],Particles[j]));
            }
        }
    }
    Interactons.Resize(0);
    for (size_t i=0; i<CInteractons.Size(); i++)
    {
        if(CInteractons[i]->UpdateContacts(Alpha)) Interactons.Push(CInteractons[i]);
    }
    for (size_t i=0; i<BInteractons.Size(); i++)
    {
        if(BInteractons[i]->UpdateContacts(Alpha)) Interactons.Push(BInteractons[i]);
    }
#endif
}

inline void Domain::ResetBoundaries()
{
    /*
    ParXmax.Resize(0);
    for (size_t i=0;i<Particles.Size();i++)
    {
        Particle * Pa = Particles[i];
        if ((Pa->MinX()>Xmax)&&Pa->IsFree())
        {
            Vec3_t v(Xmin-Xmax,0.0,0.0);
            Pa->Translate(v);
        }
        if ((Pa->MaxX()<Xmin)&&Pa->IsFree())
        {
            Vec3_t v(Xmax-Xmin,0.0,0.0);
            Pa->Translate(v);
        }
        if ((Pa->MaxX()>Xmax-2.0*Alpha-2.0*MaxDmax)&&Pa->IsFree())
        {
            ParXmax.Push(Pa);
        }
    }
    
    for (size_t i=0;i<ParXmax.Size();i++) 
    {
        Particle * P1 = ParXmax[i];
        Vec3_t v(Xmin-Xmax,0.0,0.0);
        P1->Translate(v);
        for (size_t j=0; j<Particles.Size(); j++)
        {
            Particle * P2 = Particles[j];
            if (P1==P2||!P2->IsFree()||ParXmax.Has(P2)) continue;
            bool close = (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*Alpha);
            if (!close) continue;
            std::set<std::pair<Particle *, Particle *> >::iterator it = PListofpairs.find(std::make_pair(P1,P2));
            if (it != PListofpairs.end())
            {
                continue;
            }
            PListofpairs.insert(std::make_pair(P1,P2));
            // if both particles are spheres (just one vertex)
            if (P1->Verts.Size()==1 && P2->Verts.Size()==1)
            {
                CPInteractons.Push (new CInteractonSphere(P1,P2));
            }

            // normal particles
            else
            {
                CPInteractons.Push (new CInteracton(P1,P2));
            }
        }
    }
    PInteractons.Resize(0);
    for (size_t i=0; i<CPInteractons.Size(); i++)
    {
        if(CPInteractons[i]->UpdateContacts(Alpha)) PInteractons.Push(CPInteractons[i]);
    }
    Vec3_t vmin(Xmax-Xmin,0.0,0.0);
    for (size_t i=0;i<ParXmax.Size();i++) ParXmax[i]->Translate(vmin);
    */
}

inline void Domain::CalcForceSphere()
{
    //std::cout << "Pairs size = " << ListPosPairs.Size() << std::endl;
    //std::set<std::pair<Particle *,Particle *> >::iterator it;
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
#endif
    //for (it=Listofpairs.begin();it!=Listofpairs.end();++it)
    for (size_t np=0;np<ListPosPairs.size();np++)
    {
        //std::cout << "1" << std::endl;
        size_t i = ListPosPairs[np].first;
        size_t j = ListPosPairs[np].second;
        //std::cout << i << " " << j << std::endl;
        DEM::Particle * P1 = Particles[i];
        DEM::Particle * P2 = Particles[j];
        //DEM::Particle * P1 = it->first;
        //DEM::Particle * P2 = it->second;
        //size_t i = P1->Index;
        //size_t j = P2->Index;
        bool pi = (P1->Verts.Size()==1);
        bool pj = (P2->Verts.Size()==1);
        //bool close = (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*Alpha);
        if (!(pi&&pj)) continue;
        //std::cout << i << " " << j << std::endl;

        Vec3_t xi = P1->x;
        Vec3_t xf = P2->x;
        double dist = norm(P1->x - P2->x);
        double delta = P1->Props.R + P2->Props.R - dist;
        if (delta>0)
        {
            //std::cout << "2" << std::endl;
            Vec3_t n = (xf-xi)/dist;
            Vec3_t x = xi+n*((P1->Props.R*P1->Props.R-P2->Props.R*P2->Props.R+dist*dist)/(2*dist));
            Vec3_t t1,t2,x1,x2;
            Rotation(P1->w,P1->Q,t1);
            Rotation(P2->w,P2->Q,t2);
            x1 = x - P1->x;
            x2 = x - P2->x;
            Vec3_t vrel = -((P2->v-P1->v)+cross(t2,x2)-cross(t1,x1));
            Vec3_t vt = vrel - dot(n,vrel)*n;

            double Kn = ReducedValue(P1->Props.Kn,P2->Props.Kn);
            double Kt = 2*ReducedValue(P1->Props.Kt,P2->Props.Kt);
            double me = ReducedValue(P1->Props.m ,P2->Props.m );
            double Gn = 2*ReducedValue(P1->Props.Gn,P2->Props.Gn);
            double Gt = 2*ReducedValue(P1->Props.Gt,P2->Props.Gt);
            double beta = 2*ReducedValue(P1->Props.Beta,P2->Props.Beta);
            double eta  = 2*ReducedValue(P1->Props.Eta,P2->Props.Eta);
            //double Bn   = 2*ReducedValue(P1->Props.Bn,P2->Props.Bn);
            //double Bt   = 2*ReducedValue(P1->Props.Bt,P2->Props.Bt);
            //double eps  = 2*ReducedValue(P1->Props.eps,P2->Props.eps);
            double Mu;

            if (P1->Props.Mu>1.0e-12&&P2->Props.Mu>1.0e-12)
            {
                Mu          = std::max(P1->Props.Mu,P2->Props.Mu);
            }
            else 
            {
                Mu          = 0.0;
            }

            if (Gn < -0.001)
            {
                if (fabs(Gn)>1.0) throw new Fatal("CInteractonSphere the restitution coefficient is greater than 1");
                Gn = 2.0*sqrt((pow(log(-Gn),2.0)*(Kn/me))/(M_PI*M_PI+pow(log(-Gn),2.0)));
                Gt = 0.0;
            }
            Gn *= me;
            Gt *= me;

            Vec3_t Fn = Kn*delta*n;
            //std::pair<int,int> p;
            //p = std::make_pair(i,j);
            size_t p = HashFunction(i,j);
            if (FricSpheres.count(p)==0) 
            {             
#ifdef USE_OMP
                omp_set_lock  (&lck);
#endif
                FricSpheres[p] = OrthoSys::O;
#ifdef USE_OMP
                omp_unset_lock(&lck);
#endif
            }
            FricSpheres[p] += vt*Dt;
            FricSpheres[p] -= dot(FricSpheres[p],n)*n;
            Vec3_t tan = FricSpheres[p];
            if (norm(tan)>0.0) tan/=norm(tan);
            if (norm(FricSpheres[p])>Mu*norm(Fn)/Kt)
            {
                FricSpheres[p] = Mu*norm(Fn)/Kt*tan;
            }
            Vec3_t F = Fn + Kt*FricSpheres[p] + Gn*dot(n,vrel)*n + Gt*vt;
            //Vec3_t F = Fn + P1->Props.Gn*dot(n,vrel)*n + P1->Props.Gt*vt;
            Vec3_t F1   = -F;
            Vec3_t F2   =  F;

            Vec3_t T, Tt;
            Tt = cross (x1,F);
            Quaternion_t q;
            Conjugate (P1->Q,q);
            Rotation  (Tt,q,T);
            Vec3_t T1 = -T;

            Tt = cross (x2,F);
            Conjugate (P2->Q,q);
            Rotation  (Tt,q,T);
            Vec3_t T2 =  T;
            
            //std::cout << "3" << std::endl;


            //rolling resistance
            if (dot(F,n)<0) F-=dot(F,n)*n;
            Vec3_t Normal = Fn/norm(Fn);
            Vec3_t Vr = P1->Props.R*P2->Props.R*cross(Vec3_t(t1 - t2),Normal)/(P1->Props.R+P2->Props.R);
            if (RollSpheres.count(p)==0) 
            {
#ifdef USE_OMP
                omp_set_lock  (&lck);
#endif
                RollSpheres[p] = OrthoSys::O;
#ifdef USE_OMP
                omp_unset_lock(&lck);
#endif
            }
            RollSpheres[p] += Vr*Dt;
            RollSpheres[p] -= dot(RollSpheres[p],Normal)*Normal;
            tan = RollSpheres[p];
            if (norm(tan)>0.0) tan/=norm(tan);
            double Kr = beta*Kt;
            if (norm(RollSpheres[p])>eta*Mu*norm(Fn)/Kr)
            {
                RollSpheres[p] = eta*Mu*norm(Fn)/Kr*tan;
            }
            Vec3_t Ft = -Kr*RollSpheres[p];

            //std::cout << "4" << std::endl;

            Tt = P1->Props.R*cross(Normal,Ft);
            Conjugate (P1->Q,q);
            Rotation  (Tt,q,T);
            T1 += T;
//
            Tt = P2->Props.R*cross(Normal,Ft);
            Conjugate (P2->Q,q);
            Rotation  (Tt,q,T);
            T2 -= T;

#ifdef USE_OMP
            omp_set_lock  (&P1->lck);
#endif
            P1->F += F1;
            P1->T += T1;
#ifdef USE_OMP
            omp_unset_lock(&P1->lck);
            omp_set_lock  (&P2->lck);
#endif
            P2->F += F2;
            P2->T += T2;
#ifdef USE_OMP
            omp_unset_lock(&P2->lck);
#endif
        }
    }
    
}

// Auxiliar methods

inline void Domain::LinearMomentum (Vec3_t & L)
{
    L = 0.,0.,0.;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        L += Particles[i]->Props.m*Particles[i]->v;
    }
}

inline void Domain::AngularMomentum (Vec3_t & L)
{
    L = 0.,0.,0.;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Vec3_t t1,t2;
        t1 = Particles[i]->I(0)*Particles[i]->w(0),Particles[i]->I(1)*Particles[i]->w(1),Particles[i]->I(2)*Particles[i]->w(2);
        Rotation (t1,Particles[i]->Q,t2);
        L += Particles[i]->Props.m*cross(Particles[i]->x,Particles[i]->v)+t2;
    }
}

inline double Domain::CalcEnergy (double & Ekin, double & Epot)
{
    // kinematic energy
    Ekin = 0.0;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Ekin += 0.5*Particles[i]->Props.m*dot(Particles[i]->v,Particles[i]->v)
                + 0.5*(Particles[i]->I(0)*Particles[i]->w(0)*Particles[i]->w(0)
                      +Particles[i]->I(1)*Particles[i]->w(1)*Particles[i]->w(1)
                      +Particles[i]->I(2)*Particles[i]->w(2)*Particles[i]->w(2));
    }

    // potential energy
    Epot = 0.0;
    for (size_t i=0; i<CInteractons.Size(); i++)
    {
        Epot += CInteractons[i]->Epot;
    }

    // total energy
    return Ekin + Epot;
}

inline double Domain::CriticalDt ()
{
    double MaxKn   =  0.0;
    double MaxBn   =  0.0;
    double MinMass = -1.0;
    for (size_t i=0; i<Particles.Size(); i++) 
    { 
        if (Particles[i]->IsFree())
        {
            if (Particles[i]->Props.Kn > MaxKn  ) MaxKn   = Particles[i]->Props.Kn;
            if (Particles[i]->Props.m  < MinMass||(MinMass<0.0)) MinMass = Particles[i]->Props.m;
        }
    }
    for (size_t i=0; i<BInteractons.Size(); i++)
    {
        double pbn = std::max(BInteractons[i]->Bn/BInteractons[i]->L0,BInteractons[i]->Bt/BInteractons[i]->L0);
        if (pbn > MaxBn) MaxBn = pbn;
    }

    return 0.1*sqrt(MinMass/(MaxKn+MaxBn));
}

inline void Domain::EnergyOutput (size_t IdxOut, std::ostream & OF)
{
    // header
    if (IdxOut==0)
    {
        OF << Util::_10_6 << "Time" << Util::_8s << "Ekin" << Util::_8s << "Epot" << Util::_8s << "Evis" << Util::_8s << "Efric" << Util::_8s << "Wext" << std::endl;
    }
    double Ekin,Epot;
    CalcEnergy(Ekin,Epot);
    OF << Util::_10_6 << Time << Util::_8s << Ekin << Util::_8s << Epot << Util::_8s << Evis << Util::_8s << Efric << Util::_8s << Wext << std::endl;
}

inline void Domain::GetGSD (Array<double> & X, Array<double> & Y, Array<double> & D, size_t NDiv) const
{
    // calc GSD information
    Array<double> Vg;
    double Vs = 0.0;

    for (size_t i=0; i<Particles.Size(); i++)
    {
        Particle * P = Particles[i];
        double Diam = sqrt((P->MaxX()-P->MinX())*(P->MaxX()-P->MinX())+(P->MaxY()-P->MinY())*(P->MaxY()-P->MinY())+(P->MaxZ()-P->MinZ())*(P->MaxX()-P->MinX()));
        Vs += Particles[i]->Props.V;
        Vg.Push(Particles[i]->Props.V);
        D.Push(Diam);
    }
    double Dmin  = D[D.TheMin()]; // minimum diameter
    double Dmax  = D[D.TheMax()]; // maximum diameter
    double Dspan = (Dmax-Dmin)/NDiv;
    for (size_t i=0; i<=NDiv; i++)
    {
        X.Push (i*Dspan+Dmin);
        double cumsum = 0;
        for (size_t j=0; j<D.Size(); j++)
        {
            if (D[j]<=i*Dspan+Dmin) cumsum++;
        }
        Y.Push (cumsum/Particles.Size());
    }
}

inline void Domain::Clusters ()
{
    Array<int> connections;
    for (size_t i=0;i<BInteractons.Size();i++)
    {
        if (BInteractons[i]->valid)
        {
            connections.Push(BInteractons[i]->P1->Index);
            connections.Push(BInteractons[i]->P2->Index);
        }
    }

    Util::Tree tree(connections);
    tree.GetClusters(Listofclusters);
    for (size_t i=0;i<Listofclusters.Size();i++)
    {
        for (size_t j=0;j<Listofclusters[i].Size();j++)
        {
            Particles[Listofclusters[i][j]]->Cluster = i;
        }
    }
    //std::cout << Listofclusters.Size() << std::endl;
    //std::cout << BInteractons.Size() << std::endl;
}

inline void Domain::DelParticles (Array<int> const & Tags)
{
    Array<int> idxs; // indices to be deleted
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        for (size_t j=0; j<Tags.Size(); ++j)
        {
            if (Particles[i]->Tag==Tags[j]) idxs.Push(i);
        }
    }
    if (idxs.Size()<1) throw new Fatal("Domain::DelParticles: Could not find any particle to be deleted");
    Particles.DelItems (idxs);
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        Particles[i]->Index = i;
    }
}

inline Particle * Domain::GetParticle (int Tag, bool Check)
{
    size_t idx   = 0;
    size_t count = 0;
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        if (Particles[i]->Tag==Tag)
        {
            if (!Check) return Particles[i];
            idx = i;
            count++;
        }
    }
    if      (count==0) throw new Fatal("Domain::GetParticle: Could not find Particle with Tag==%d",Tag);
    else if (count>1)  throw new Fatal("Domain::GetParticle: There are more than one particle with Tag==%d",Tag);
    return Particles[idx];
}

inline Particle const & Domain::GetParticle (int Tag, bool Check) const
{
    size_t idx   = 0;
    size_t count = 0;
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        if (Particles[i]->Tag==Tag)
        {
            if (!Check) return (*Particles[i]);
            idx = i;
            count++;
        }
    }
    if      (count==0) throw new Fatal("Domain::GetParticle: Could not find Particle with Tag==%d",Tag);
    else if (count>1)  throw new Fatal("Domain::GetParticle: There are more than one particle with Tag==%d",Tag);
    return (*Particles[idx]);
}

inline void Domain::GetParticles (int Tag, Array<Particle*> & P)
{
    P.Resize(0);
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        if (Particles[i]->Tag==Tag) P.Push(Particles[i]);
    }
    if (P.Size()==0) throw new Fatal("Domain::GetParticles: Could not find any Particle with Tag==%d",Tag);
}

}; // namespace DEM

#endif // MECHSYS_DEM_DOMAIN_H
