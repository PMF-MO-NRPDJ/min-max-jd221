#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <limits>

#include "dune/common/parallel/mpihelper.hh"
#include <dune/common/exceptions.hh> 

#include <dune/grid/uggrid.hh>  
#include <dune/grid/common/gridinfo.hh> 
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

// Izračunaj kut između p1-p2 i p3-p2
template <typename Point>
double kut(Point const & p1, Point const & p2, Point const & p3)
{
   // računanje kuta (VAŠ KOD) ide ovdje    
   // Point će biti Dune::FieldVector<double, dim>. Norma se računa 
   // pomoću funkcije članice klase two_norm(), a skalarni produkt 
   // pomoću funkcije članice klase dot(). Pogledati dokumentaciju klase
   //  Dune::FieldVector<double, dim>.
    return kut;
}

int main(int argc, char** argv)
{
    const int dim = 2;
    using GridType = Dune::UGGrid<dim>;
    using LeafGridView = ...;
     
    // UČITATI 2D MREŽU IZ GMSH DATOTEKE KOJA JE ZADANA KAO ARGUMENT KOMANDNE LINIJE.
    
    // gridView JE LeafGridView!
    for(auto const & element : elements(gridView))
    {

/*     VAŠ KOD dolazi ovdje.
 *     RAČUNATI MIN I MAX KUT U SVAKOM ELEMENTU. 
*/
    } 

   // ISPISATI BROJ ELEMENATA; MINIMALNI I MAKSIMALNI KUT U STUPNJEVIMA:
    

    // Ispis mreže u VTK formatu (u datoteku poluvijenac.vtu)
    Dune::VTKWriter<LeafGridView> vtkwriter(gridView);
    vtkwriter.write("poluvijenac");

    return 0;
}
