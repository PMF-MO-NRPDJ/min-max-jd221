#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <vector>
#include <cmath>


#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>


// Funkcija racuna kut T1T2T3
template <typename Point>
double angle(Point const & T1, Point const & T2, Point const & T3)
{
    Dune::FieldVector<double,2> vector1=T1-T2;
    Dune::FieldVector<double,2> vector2=T3-T2;

    // cos^-1(x*y/(norma(x)*norma(y)))

    // FieldVector funkcije
    double skalarni_produkt=dot(vector1,vector2);
    double norme=vector1.two_norm()*vector2.two_norm();

    double kut=std::acos(skalarni_produkt/norme);

    // Stupnjevi
    kut=180*kut/M_PI; // M_PI=pi

    return kut;
}

int main(int argc, char** argv)
{
    if (argc == 1) {
    std::cerr << "Usage: " << argv[0] << " ime_grid_datoteke.msh"
                  << std::endl;
    std::exit(1);
    }

    Dune::MPIHelper::instance(argc, argv);

    const int dim = 2;
    using GridType = Dune::UGGrid<dim>;
    using GridView = GridType::LeafGridView;

    bool verbosity = true;
    bool insertBoundarySegments = false;

    std::unique_ptr<GridType> p_grid = Dune::GmshReader<GridType>::read(argv[1], verbosity, insertBoundarySegments);
    auto gridView = p_grid->leafGridView();

    using GlobCoo = GridView::template Codim<0>::Geometry::GlobalCoordinate;


    double kut_min=400;
    double kut_max=-400;
    int index_min=0;
    int index_max=0;

    int count=0;
    double kut;

    for(auto const & element : elements(gridView))
    {
        // Geometrijski tip elementa (tip Dune::GeometryType)
        auto gt = element.type();
        auto geom=element.geometry();
        // Broj vrhova elementa
        auto n_v = geom.corners(); // trokutovi --> =3

        // Uzmi koordinate svih vrhova elemeta
        std::vector<GlobCoo> coo_v(n_v);

        kut=angle(element.geometry().corner(0), element.geometry().corner(1), element.geometry().corner(2));

        if(kut<kut_min){
          kut_min=kut;
          index_min=count;
        }

        if(kut>kut_max){
          kut_max=kut;
          index_max=count;
        }

        kut=angle(element.geometry().corner(1), element.geometry().corner(2), element.geometry().corner(0));

        if(kut<kut_min){
          kut_min=kut;
          index_min=count;
        }

        if(kut>kut_max){
          kut_max=kut;
          index_max=count;
        }

        kut=angle(element.geometry().corner(2), element.geometry().corner(0), element.geometry().corner(1));

        if(kut<kut_min){
          kut_min=kut;
          index_min=count;
        }

        if(kut>kut_max){
          kut_max=kut;
          index_max=count;
        }

        ++count;
    }


    //Rješenje - min i max kut:
    std::cout<<"Ukupni broj elemenata: "<<count<<"\n";
    std::cout<<"Najmanji kut: "<<kut_min<<" -> element "<<index_min<<"\n";
    std::cout<<"Najveci kut: "<<kut_max<<" -> element "<<index_max<<"\n";

    std::cout<<std::endl;

    Dune::VTKWriter<GridView> vtkwriter(gridView);
    vtkwriter.write("poluvijenac");

    return 0;
}
