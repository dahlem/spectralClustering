#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "arerror.h"
#include "arlsmat.h"
#include "arlgsym.h"
#include "lsymsol.h"

int main(int argc, char* argv[])
{
    int nev = boost::lexical_cast<int>(argv[6]);
    int ncv = 0;
    int maxit = 0;
    int nconv = 0;
    double tol = 0.00001;
    char* which = "SM";
    char up = 'U';

    std::string eigenValues = std::string(argv[1]) + "/eigenValues.csv";
    std::string eigenVectors = std::string(argv[1]) + "/eigenVectors.csv";

    // Compressed Row Storage (CRS) storage constraints:
    // 2 * nnz + n + 1
    // D has as many non-zeros as rows (or columns)
    int nnzD = boost::lexical_cast<int>(argv[3]); // non-zeros in D
    double *nzvalD = new double[nnzD];
    int *irowD = new int[nnzD];
    int *pcolD = new int[nnzD+1];

    int nnzL = boost::lexical_cast<int>(argv[2]) + nnzD; // non-zeros in L
    double *nzvalL = new double[nnzL];
    int *irowL = new int[nnzL];
    int *pcolL = new int[nnzD+1];

    std::ifstream DFile;
    DFile.open(argv[4]);
    if (!DFile.is_open()) {
        std::cerr << "Could not open file: " << argv[4] << std::endl;
        return EXIT_FAILURE;
    }

    std::ifstream WFile;
    WFile.open(argv[5]);
    if (!WFile.is_open()) {
        std::cerr << "Could not open file: " << argv[5] << std::endl;
        return EXIT_FAILURE;
    }
    
    // parse
    std::string line;
    // read diagonal matrix
    while (!DFile.eof()) {
        std::getline(DFile, line);
        boost::trim(line);

        if (line != "") {
            std::vector<std::string> splitVec;
            boost::split(splitVec, line, boost::is_any_of(","), boost::token_compress_on);
            boost::uint32_t d = boost::lexical_cast<boost::uint32_t>(splitVec[0]);
            double v = boost::lexical_cast<double>(splitVec[1]);
            nzvalD[d-1] = v;
            irowD[d-1] = d-1;
            pcolD[d-1] = d-1;
        }
    }
    pcolD[nnzD] = nnzD;

    // read W matrix and store into a laplacian matrix L=D-W
    boost::uint32_t idx = 0;
    boost::uint32_t curCol = 0;
    pcolL[curCol] = idx;
    
    while (!WFile.eof()) {
        std::getline(WFile, line);
        boost::trim(line);

        if (line != "") {
            std::vector<std::string> splitVec;
            boost::split(splitVec, line, boost::is_any_of(","), boost::token_compress_on);
            boost::uint32_t r = boost::lexical_cast<boost::uint32_t>(splitVec[0]);
            boost::uint32_t c = boost::lexical_cast<boost::uint32_t>(splitVec[1]);
            double v = boost::lexical_cast<double>(splitVec[2]);
            if ((c-1) != curCol) {
                // fill in the diagonal(s) of previous column(s)
                while (curCol < (c-1)) {
                    nzvalL[idx] = nzvalD[curCol];
                    irowL[idx] = curCol;

                    idx++;
                    curCol++; // start new column
                    pcolL[curCol] = idx; // index of new column into nzvalL and irowL
                }
            }
            // we are operating on the same column as before
            nzvalL[idx] = -v;
            irowL[idx] = r-1;
            idx++;
        }
    }
    // fill the necessary diagonal elements
    while (curCol < nnzD) {
        nzvalL[idx] = nzvalD[curCol];
        irowL[idx] = curCol;

        idx++;
        curCol++; // start new column
        pcolL[curCol] = idx; // index of new column into nzvalL and irowL
    }
    
    pcolL[nnzD] = idx;
    WFile.close();
    
    ARluSymMatrix<double> D(nnzD, nnzD, nzvalD, irowD, pcolD, up);
    ARluSymMatrix<double> L(nnzD, nnzL, nzvalL, irowL, pcolL, up);
    ARluSymGenEig<double> dprob(nev, L, D, which, ncv, tol, maxit);

    dprob.Trace();
    std::cout << "Compute eigensystem..." << std::endl;
    try {
        nconv = dprob.FindEigenvectors();
    } catch (ArpackError) {
        return 0;
    }

    if (dprob.EigenvaluesFound()) {
        std::ofstream outVals(eigenValues.c_str(), std::ios::out);
        for (int i = 0; i < nconv; ++i) {
            outVals << dprob.Eigenvalue(i) << std::endl;
        }
        outVals.close();

        // serialise the eigenvectors for k-means clustering
        std::ofstream outVecs(eigenVectors.c_str(), std::ios::out);
        for (int j = 0; j < nnzD; ++j) { 
            for (int i = 0; i < (nconv-1); ++i) { // let's output all vectors
                outVecs << dprob.Eigenvector(i, j) << ",";
            }
            outVecs << dprob.Eigenvector(nconv-1,j) << std::endl;
        }
        outVecs.close();
    }
    
    Solution(L, D, dprob);

    delete[] nzvalD;
    delete[] irowD;
    delete[] pcolD;
    delete[] nzvalL;
    delete[] irowL;
    delete[] pcolL;
}
