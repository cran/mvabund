// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// Rcpp.h: R/C++ interface class library
//
// Copyright (C) 2005 - 2006 Dominick Samperi
// Copyright (C) 2008 - 2009 Dirk Eddelbuettel
//
// This library is free software; you can redistribute it and/or modify it 
// under the terms of the GNU Lesser General Public License as published by 
// the Free Software Foundation; either version 2.1 of the License, or (at 
// your option) any later version.
//
// This library is distributed in the hope that it will be useful, but 
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public 
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public License 
// along with this library; if not, write to the Free Software Foundation, 
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 

#ifndef mvRcpp_hpp
#define mvRcpp_hpp

#include <iostream>
#include <sstream>
#include <string>
#include <list>
#include <map>
#include <stdexcept>
#include <vector>

// include R headers, but set R_NO_REMAP and access everything via Rf_ prefixes
#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>

// #ifdef BUILDING_DLL
// #define RcppExport extern "C" __declspec(dllexport)
// #else
#define RcppExport extern "C"
// #endif

char *copyMessageToR(const char* const mesg);

class RcppParams {
public:
    RcppParams(SEXP params);
    // functions
    void   	checkNames(char *inputNames[], int len);
    bool	exists(std::string name);
    double 	getDoubleValue(std::string name);
    int    	getIntValue(std::string name);
    bool   	getBoolValue(std::string name);
    std::string getStringValue(std::string name);
private:
    std::map<std::string, int> pmap;
    SEXP _params;
};

class RcppNumList {
public:
    RcppNumList(SEXP theList) {
	if (!Rf_isNewList(theList))
	    throw std::range_error("RcppNumList: non-list passed to constructor");
        len = Rf_length(theList);
        names = Rf_getAttrib(theList, R_NamesSymbol);
        namedList = theList;
    }
    std::string getName(int i) {
        if (i < 0 || i >= len) {
	    std::ostringstream oss;
	    oss << "RcppNumList::getName: index out of bounds: " << i;
	    throw std::range_error(oss.str());
	}
        return std::string(CHAR(STRING_ELT(names,i)));
    }
    double getValue(int i) {
        if (i < 0 || i >= len) {
	    std::ostringstream oss;
	    oss << "RcppNumList::getValue: index out of bounds: " << i;
	    throw std::range_error(oss.str());
	}
	SEXP elt = VECTOR_ELT(namedList, i);
	if (Rf_isReal(elt))
	    return REAL(elt)[0];
	else if (Rf_isInteger(elt))
	    return (double)INTEGER(elt)[0];
	else
	    throw std::range_error("RcppNumList: contains non-numeric value");
	return 0; // never get here
    }
    int size() { return len; }
private:
    int len;
    SEXP namedList;
    SEXP names;
};

template <typename T>
class RcppVector {
public:
    RcppVector(SEXP vec);
    RcppVector(int len);
    int size() { return len; }
    inline T& operator()(int i) {
	if (i < 0 || i >= len) {
	    std::ostringstream oss;
	    oss << "RcppVector: subscript out of range: " << i;
	    throw std::range_error(oss.str());
	}
	return v[i];
    }
    T *cVector();
    std::vector<T> stlVector();
private:
    int len;
    T *v;
};

template <typename T>
class RcppMatrix {
public:
    RcppMatrix(SEXP mat);
    RcppMatrix(int nx, int ny);
    int getDim1() { return dim1; }
    int getDim2() { return dim2; }
    inline T& operator()(int i, int j) {
	if (i < 0 || i >= dim1 || j < 0 || j >= dim2) {
	    std::ostringstream oss;
	    oss << "RcppMatrix: subscripts out of range: " << i << ", " << j;
	    throw std::range_error(oss.str());
	}
	return a[i][j];
    }
    T **cMatrix();
    std::vector<std::vector<T> > stlMatrix();
private:
    int dim1, dim2;
    T **a;
};


template <typename T>
class RcppVectorView {
public:
    RcppVectorView(SEXP vec);
    inline int size() const { return len; }
    inline T operator()(int i) const {
	if (i < 0 || i >= len) {
	    std::ostringstream oss;
	    oss << "RcppVectorView: subscript out of range: " << i;
	    throw std::range_error(oss.str());
	}
	return v[i];
    }
private:
    int len;
    T *v;
};

template <typename T>
class RcppMatrixView {
public:
    RcppMatrixView(SEXP mat);
    inline int dim1() const { return d1; }
    inline int dim2() const { return d2; }
    inline T operator()(int i, int j) const {
	if (i < 0 || i >= d1 || j < 0 || j >= d2) {
	    std::ostringstream oss;
	    oss << "RcppMatrixView: subscripts out of range: " << i << ", " << j;
	    throw std::range_error(oss.str());
	}
	return a[i + d1 * j];
    }
private:
    T *a;
    int d1, d2;
};

class RcppFunction {
public:
    RcppFunction(SEXP fn) : fn(fn) { 
	if (!Rf_isFunction(fn))
	    throw std::range_error("RcppFunction: non-function where function expected");
	numProtected = 0;
	currListPosn = 0;
        listSize = 0;
        vectorArg = listArg = R_NilValue;
    };
    ~RcppFunction() {
	UNPROTECT(numProtected);
    }
    SEXP listCall();
    SEXP vectorCall();
    void setRVector(std::vector<double>& v);
    void setRListSize(int size);
    void appendToRList(std::string name, double value);
    void appendToRList(std::string name, int value);
    void clearProtectionStack() {
	UNPROTECT(numProtected);
	numProtected = 0;
    }

private:
    SEXP fn, listArg, vectorArg;
    int listSize, currListPosn, numProtected;
    std::vector<std::string> names;
};

class RcppResultSet {
public:
    RcppResultSet() { numProtected = 0; }
    void add(std::string, double);
    void add(std::string, int);
    void add(std::string, double *, int);
    void add(std::string, int *, int);
    void add(std::string, double **, int, int);
    void add(std::string, int **, int, int);
    void add(std::string, std::vector<double>&);
    void add(std::string, std::vector<int>&);
    void add(std::string, std::vector<std::vector<double> >&);
    void add(std::string, std::vector<std::vector<int> >&);
    void add(std::string, RcppVector<int>&);
    void add(std::string, RcppVector<double>&);
    void add(std::string, RcppMatrix<int>&);
    void add(std::string, RcppMatrix<double>&);
    void add(std::string, SEXP, bool isProtected);
    SEXP getReturnList();
protected:
    int numProtected;
    std::list<std::pair<std::string,SEXP> > values;
};


#endif
