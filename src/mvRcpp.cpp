// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// Rcpp.cpp: R/C++ interface class library
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

#include "mvRcpp.h"

RcppParams::RcppParams(SEXP params) {
    if (!Rf_isNewList(params))
	throw std::range_error("RcppParams: non-list passed to constructor");
    int len = Rf_length(params);
    SEXP names = Rf_getAttrib(params, R_NamesSymbol);
    if (names == R_NilValue)
	throw std::range_error("RcppParams: list must have named elements");
    for (int i = 0; i < len; i++) {
	std::string nm = std::string(CHAR(STRING_ELT(names,i)));
	if (nm.size() == 0)
	    throw std::range_error("RcppParams: all list elements must be named");
	pmap[nm] = i;
    }
    _params = params;
}

void RcppParams::checkNames(char *inputNames[], int len) {
    for (int i = 0; i < len; i++) {
	std::map<std::string,int>::iterator iter = pmap.find(inputNames[i]);
	if (iter == pmap.end()) {
	    std::string mesg = "RcppParams::checkNames: missing required parameter ";
	    throw std::range_error(mesg+inputNames[i]);
	}
    }
}

bool RcppParams::exists(std::string name) {
    bool rc = true;
    std::map<std::string,int>::iterator iter = pmap.find(name);
    if (iter == pmap.end()) {
	rc = false;
    }
    return rc;
}

double RcppParams::getDoubleValue(std::string name) {
    std::map<std::string,int>::iterator iter = pmap.find(name);
    if (iter == pmap.end()) {
	std::string mesg = "RcppParams::getDoubleValue: no such name: ";
	throw std::range_error(mesg+name);
    }
    int posn = iter->second;
    SEXP elt = VECTOR_ELT(_params,posn);
    if (!Rf_isNumeric(elt) || Rf_length(elt) != 1) {
	std::string mesg = "RcppParams::getDoubleValue: must be scalar ";
	throw std::range_error(mesg+name);
    }
    if (Rf_isInteger(elt))
	return (double)INTEGER(elt)[0];
    else if (Rf_isReal(elt))
	return REAL(elt)[0];
    else {
	std::string mesg = "RcppParams::getDoubleValue: invalid value for ";
	throw std::range_error(mesg+name);
    }
    return 0; // never get here
}

int RcppParams::getIntValue(std::string name) {
    std::map<std::string,int>::iterator iter = pmap.find(name);
    if (iter == pmap.end()) {
	std::string mesg = "RcppParams::getIntValue: no such name: ";
	throw std::range_error(mesg+name);
    }
    int posn = iter->second;
    SEXP elt = VECTOR_ELT(_params,posn);
    if (!Rf_isNumeric(elt) || Rf_length(elt) != 1) {
	std::string mesg = "RcppParams::getIntValue: must be scalar: ";
	throw std::range_error(mesg+name);
    }
    if (Rf_isInteger(elt))
	return INTEGER(elt)[0];
    else if (Rf_isReal(elt))
	return (int)REAL(elt)[0];
    else {
	std::string mesg = "RcppParams::getIntValue: invalid value for: ";
	throw std::range_error(mesg+name);
    }
    return 0; // never get here
}

bool RcppParams::getBoolValue(std::string name) {
    std::map<std::string,int>::iterator iter = pmap.find(name);
    if (iter == pmap.end()) {
	std::string mesg = "RcppParams::getBoolValue: no such name: ";
	throw std::range_error(mesg+name);
    }
    int posn = iter->second;
    SEXP elt = VECTOR_ELT(_params,posn);
    if (Rf_isLogical(elt))
	return INTEGER(elt)[0];
    else {
	std::string mesg = "RcppParams::getBoolValue: invalid value for: ";
	throw std::range_error(mesg+name);
    }
    return false; // never get here
}

std::string RcppParams::getStringValue(std::string name) {
    std::map<std::string,int>::iterator iter = pmap.find(name);
    if (iter == pmap.end()) {
       std::string mesg = "RcppParams::getStringValue: no such name: ";
    throw std::range_error(mesg+name);
    }
    int posn = iter->second;
    SEXP elt = VECTOR_ELT(_params,posn);
    if (Rf_isString(elt))
       return std::string(CHAR(STRING_ELT(elt,0)));
    else {
       std::string mesg = "RcppParams::getStringValue: invalid value for: ";
       throw std::range_error(mesg+name);
    }

    return ""; // never get here
}

template <typename T>
RcppVector<T>::RcppVector(SEXP vec) {
    int i;

    // The function Rf_isVector returns TRUE for vectors AND
    // matrices, so it does not distinguish. We could
    // check the dim attribute here to be sure that it
    // is not present (i.e., dimAttr == R_NilValue, not 0!).
    // But it is easier to simply check if it is set via
    // Rf_isMatrix (in which case we don't have a vector).
    if (!Rf_isNumeric(vec) || Rf_isMatrix(vec) || Rf_isLogical(vec))
	throw std::range_error("RcppVector: invalid numeric vector in constructor");
    len = Rf_length(vec);
    v = (T *)R_alloc(len, sizeof(T));
    if (Rf_isInteger(vec)) {
	for (i = 0; i < len; i++)
	    v[i] = (T)(INTEGER(vec)[i]);
    }	
    else if (Rf_isReal(vec)) {
	for (i = 0; i < len; i++)
	    v[i] = (T)(REAL(vec)[i]);
    }
}

template <typename T>
RcppVector<T>::RcppVector(int _len) {
    len = _len;
    v = (T *)R_alloc(len, sizeof(T));
    for (int i = 0; i < len; i++)
	v[i] = 0;
}

template <typename T>
T *RcppVector<T>::cVector() {
    T* tmp = (T *)R_alloc(len, sizeof(T));
    for (int i = 0; i < len; i++)
	tmp[i] = v[i];
    return tmp;
}

template <typename T>
std::vector<T> RcppVector<T>::stlVector() {
    std::vector<T> tmp(len);
    for (int i = 0; i < len; i++)
	tmp[i] = v[i];
    return tmp;
}

template <typename T>
RcppMatrix<T>::RcppMatrix(SEXP mat) {

    if (!Rf_isNumeric(mat) || !Rf_isMatrix(mat))
	throw std::range_error("RcppMatrix: invalid numeric matrix in constructor");

    // Get matrix dimensions
    SEXP dimAttr = Rf_getAttrib(mat, R_DimSymbol);
    dim1 = INTEGER(dimAttr)[0];
    dim2 = INTEGER(dimAttr)[1];

    // We guard against  the possibility that R might pass an integer matrix.
    // Can be prevented using R code: temp <- as.double(a), dim(temp) <- dim(a)
    int i,j;
    int isInt = Rf_isInteger(mat);
    T *m = (T *)R_alloc(dim1*dim2, sizeof(T));
    a = (T **)R_alloc(dim1, sizeof(T *));
    for (i = 0; i < dim1; i++)
	a[i] = m + i*dim2;
    if (isInt) {
	for (i=0; i < dim1; i++)
	    for (j=0; j < dim2; j++)
		a[i][j] = (T)(INTEGER(mat)[i+dim1*j]);
    }	
    else {
	for (i=0; i < dim1; i++)
	    for (j=0; j < dim2; j++)
		a[i][j] = (T)(REAL(mat)[i+dim1*j]);
    }	
}

template <typename T>
RcppMatrix<T>::RcppMatrix(int _dim1, int _dim2) {
    dim1 = _dim1;
    dim2 = _dim2;
    int i,j;
    T *m = (T *)R_alloc(dim1*dim2, sizeof(T));
    a = (T **)R_alloc(dim1, sizeof(T *));
    for (i = 0; i < dim1; i++)
	a[i] = m + i*dim2;
    for (i=0; i < dim1; i++)
	for (j=0; j < dim2; j++)
	    a[i][j] = 0;
}

template <typename T>
std::vector<std::vector<T> > RcppMatrix<T>::stlMatrix() {
    int i,j;
    std::vector<std::vector<T> > temp;
    for (i = 0; i < dim1; i++) {
	temp.push_back(std::vector<T>(dim2));
    }
    for (i = 0; i < dim1; i++)
	for (j = 0; j < dim2; j++)
	    temp[i][j] = a[i][j];
    return temp;
}

template <typename T>
T **RcppMatrix<T>::cMatrix() {
    int i,j;
    T *m = (T *)R_alloc(dim1*dim2, sizeof(T));
    T **tmp = (T **)R_alloc(dim1, sizeof(T *));
    for (i = 0; i < dim1; i++)
	tmp[i] = m + i*dim2;
    for (i=0; i < dim1; i++)
	for (j=0; j < dim2; j++)
	    tmp[i][j] = a[i][j];
    return tmp;
}

// Explicit instantiation (required for external linkage)
template class RcppVector<int>;
template class RcppVector<double>;
template class RcppMatrix<int>;
template class RcppMatrix<double>;

template <typename T>
RcppVectorView<T>::RcppVectorView(SEXP vec) {
    if (!Rf_isNumeric(vec) || Rf_isMatrix(vec) || Rf_isLogical(vec))
	throw std::range_error("RcppVectorView: invalid numeric vector in constructor");
    len = Rf_length(vec);
    if (Rf_isInteger(vec)) v = (T *)(INTEGER(vec));
    else if (Rf_isReal(vec)) v = (T *)(REAL(vec));
}

template class RcppVectorView<int>;
template class RcppVectorView<double>;

template <typename T>
RcppMatrixView<T>::RcppMatrixView(SEXP mat) {
    if (!Rf_isNumeric(mat) || !Rf_isMatrix(mat))
	throw std::range_error("RcppMatrixView: invalid numeric matrix in constructor");
    // Get matrix dimensions
    SEXP dimAttr = Rf_getAttrib(mat, R_DimSymbol);
    d1 = INTEGER(dimAttr)[0];
    d2 = INTEGER(dimAttr)[1];
    if (Rf_isInteger(mat)) a = (T *)(INTEGER(mat));
    else if (Rf_isReal(mat)) a = (T *)(REAL(mat));
}

template class RcppMatrixView<int>;
template class RcppMatrixView<double>;

void RcppResultSet::add(std::string name, double x) {
    SEXP value = PROTECT(Rf_allocVector(REALSXP, 1));
    numProtected++;
    REAL(value)[0] = x;
    values.push_back(make_pair(name, value));
}

void RcppResultSet::add(std::string name, int i) {
    SEXP value = PROTECT(Rf_allocVector(INTSXP, 1));
    numProtected++;
    INTEGER(value)[0] = i;
    values.push_back(make_pair(name, value));
}

void RcppResultSet::add(std::string name, double *vec, int len) {
    if (vec == 0)
	throw std::range_error("RcppResultSet::add: NULL double vector");
    SEXP value = PROTECT(Rf_allocVector(REALSXP, len));
    numProtected++;
    for (int i = 0; i < len; i++)
	REAL(value)[i] = vec[i];
    values.push_back(make_pair(name, value));
}

void RcppResultSet::add(std::string name, int *vec, int len) {
    if (vec == 0)
	throw std::range_error("RcppResultSet::add: NULL int vector");
    SEXP value = PROTECT(Rf_allocVector(INTSXP, len));
    numProtected++;
    for (int i = 0; i < len; i++)
	INTEGER(value)[i] = vec[i];
    values.push_back(make_pair(name, value));
}

void RcppResultSet::add(std::string name, double **mat, int nx, int ny) {
    if (mat == 0)
	throw std::range_error("RcppResultSet::add: NULL double matrix");
    SEXP value = PROTECT(Rf_allocMatrix(REALSXP, nx, ny));
    numProtected++;
    for (int i = 0; i < nx; i++)
	for (int j = 0; j < ny; j++)
	    REAL(value)[i + nx*j] = mat[i][j];
    values.push_back(make_pair(name, value));
}

void RcppResultSet::add(std::string name, int **mat, int nx, int ny) {
    if (mat == 0)
	throw std::range_error("RcppResultSet::add: NULL int matrix");
    SEXP value = PROTECT(Rf_allocMatrix(INTSXP, nx, ny));
    numProtected++;
    for (int i = 0; i < nx; i++)
	for (int j = 0; j < ny; j++)
	    INTEGER(value)[i + nx*j] = mat[i][j];
    values.push_back(make_pair(name, value));
}

void RcppResultSet::add(std::string name, std::vector<int>& vec) {
    if (vec.size() == 0)
	throw std::range_error("RcppResultSet::add; zero length vector<int>");
    int len = (int)vec.size();
    SEXP value = PROTECT(Rf_allocVector(INTSXP, len));
    numProtected++;
    for (int i = 0; i < len; i++)
	INTEGER(value)[i] = vec[i];
    values.push_back(make_pair(name, value));
}

void RcppResultSet::add(std::string name, std::vector<double>& vec) {
    if (vec.size() == 0)
	throw std::range_error("RcppResultSet::add; zero length vector<double>");
    int len = (int)vec.size();
    SEXP value = PROTECT(Rf_allocVector(REALSXP, len));
    numProtected++;
    for (int i = 0; i < len; i++)
	REAL(value)[i] = vec[i];
    values.push_back(make_pair(name, value));
}

void RcppResultSet::add(std::string name, std::vector<std::vector<int> >& mat) {
    if (mat.size() == 0)
	throw std::range_error("RcppResultSet::add: zero length vector<vector<int> >");
    else if (mat[0].size() == 0)
	throw std::range_error("RcppResultSet::add: no columns in vector<vector<int> >");
    int nx = (int)mat.size();
    int ny = (int)mat[0].size();
    SEXP value = PROTECT(Rf_allocMatrix(INTSXP, nx, ny));
    numProtected++;
    for (int i = 0; i < nx; i++)
	for (int j = 0; j < ny; j++)
	    INTEGER(value)[i + nx*j] = mat[i][j];
    values.push_back(make_pair(name, value));
}

void RcppResultSet::add(std::string name, std::vector<std::vector<double> >& mat) {
    if (mat.size() == 0)
	throw std::range_error("RcppResultSet::add: zero length vector<vector<double> >");
    else if (mat[0].size() == 0)
	throw std::range_error("RcppResultSet::add: no columns in vector<vector<double> >");
    int nx = (int)mat.size();
    int ny = (int)mat[0].size();
    SEXP value = PROTECT(Rf_allocMatrix(REALSXP, nx, ny));
    numProtected++;
    for (int i = 0; i < nx; i++)
	for (int j = 0; j < ny; j++)
	    REAL(value)[i + nx*j] = mat[i][j];
    values.push_back(make_pair(name, value));
}

void RcppResultSet::add(std::string name, RcppVector<int>& vec) {
    int len = vec.size();
    int *a = vec.cVector();
    SEXP value = PROTECT(Rf_allocVector(INTSXP, len));
    numProtected++;
    for (int i = 0; i < len; i++)
	INTEGER(value)[i] = a[i];
    values.push_back(make_pair(name, value));
}

void RcppResultSet::add(std::string name, RcppVector<double>& vec) {
    int len = vec.size();
    double *a = vec.cVector();
    SEXP value = PROTECT(Rf_allocVector(REALSXP, len));
    numProtected++;
    for (int i = 0; i < len; i++)
	REAL(value)[i] = a[i];
    values.push_back(make_pair(name, value));
}

void RcppResultSet::add(std::string name, RcppMatrix<int>& mat) {
    int nx = mat.getDim1();
    int ny = mat.getDim2();
    int **a = mat.cMatrix();
    SEXP value = PROTECT(Rf_allocMatrix(INTSXP, nx, ny));
    numProtected++;
    for (int i = 0; i < nx; i++)
	for (int j = 0; j < ny; j++)
	    INTEGER(value)[i + nx*j] = a[i][j];
    values.push_back(make_pair(name, value));
}

void RcppResultSet::add(std::string name, RcppMatrix<double>& mat) {
    int nx = mat.getDim1();
    int ny = mat.getDim2();
    double **a = mat.cMatrix();
    SEXP value = PROTECT(Rf_allocMatrix(REALSXP, nx, ny));
    numProtected++;
    for (int i = 0; i < nx; i++)
	for (int j = 0; j < ny; j++)
	    REAL(value)[i + nx*j] = a[i][j];
    values.push_back(make_pair(name, value));
}


void RcppResultSet::add(std::string name, SEXP sexp, bool isProtected) {
    values.push_back(make_pair(name, sexp));
    if (isProtected)
	numProtected++;
}

SEXP RcppResultSet::getReturnList() {
    int nret = (int)values.size();
    SEXP rl = PROTECT(Rf_allocVector(VECSXP,nret));
    SEXP nm = PROTECT(Rf_allocVector(STRSXP,nret));
    std::list<std::pair<std::string,SEXP> >::iterator iter = values.begin();
    for (int i = 0; iter != values.end(); iter++, i++) {
	SET_VECTOR_ELT(rl, i, iter->second);
	SET_STRING_ELT(nm, i, Rf_mkChar(iter->first.c_str()));
    }
    Rf_setAttrib(rl, R_NamesSymbol, nm);
    UNPROTECT(numProtected+2);
    return rl;
}

SEXP RcppFunction::listCall() {
    if (names.size() != (unsigned)listSize)
	throw std::range_error("RcppFunction::listCall: no. of names != no. of items");
    if (currListPosn != listSize)
	throw std::range_error("RcppFunction::listCall: list has incorrect size");
    SEXP nm = PROTECT(Rf_allocVector(STRSXP,listSize));
    numProtected++;
    for (int i=0; i < listSize; i++)
	SET_STRING_ELT(nm, i, Rf_mkChar(names[i].c_str()));
    Rf_setAttrib(listArg, R_NamesSymbol, nm);
    SEXP R_fcall;
    PROTECT(R_fcall = Rf_lang2(fn, R_NilValue));
    numProtected++;
    SETCADR(R_fcall, listArg);
    SEXP result = Rf_eval(R_fcall, R_NilValue);
    names.clear();
    listSize = currListPosn = 0; // Ready for next call.
    return result;
}

SEXP RcppFunction::vectorCall() {
    if (vectorArg == R_NilValue)
	throw std::range_error("RcppFunction::vectorCall: vector has not been set");
    SEXP R_fcall;
    PROTECT(R_fcall = Rf_lang2(fn, R_NilValue));
    numProtected++;
    SETCADR(R_fcall, vectorArg);
    SEXP result = Rf_eval(R_fcall, R_NilValue);
    vectorArg = R_NilValue; // Ready for next call.
    return result;
}

void RcppFunction::setRVector(std::vector<double>& v) {
    vectorArg = PROTECT(Rf_allocVector(REALSXP,v.size()));
    numProtected++;
    for (int i=0; i < (int)v.size(); i++)
	REAL(vectorArg)[i] = v[i];
}

void RcppFunction::setRListSize(int n) {
    listSize = n;
    listArg = PROTECT(Rf_allocVector(VECSXP, n));
    numProtected++;
}

void RcppFunction::appendToRList(std::string name, double value) {
    if (currListPosn < 0 || currListPosn >= listSize)
	throw std::range_error("RcppFunction::appendToRList(double): list posn out of range");
    SEXP valsxp = PROTECT(Rf_allocVector(REALSXP,1));
    numProtected++;
    REAL(valsxp)[0] = value;
    SET_VECTOR_ELT(listArg, currListPosn++, valsxp);
    names.push_back(name);
}

void RcppFunction::appendToRList(std::string name, int value) {
    if (currListPosn < 0 || currListPosn >= listSize)
	throw std::range_error("RcppFunction::appendToRlist(int): posn out of range");
    SEXP valsxp = PROTECT(Rf_allocVector(INTSXP,1));
    numProtected++;
    INTEGER(valsxp)[0] = value;
    SET_VECTOR_ELT(listArg, currListPosn++, valsxp);
    names.push_back(name);
}

#include <cstring>

// Paul Roebuck has observed that the memory used by an exception message
// is not reclaimed if error() is called inside of a catch block (due to
// a setjmp() call), and he suggested the following work-around.
char *copyMessageToR(const char* const mesg) {
    char* Rmesg;
    const char* prefix = "Exception: ";
    void* Rheap = R_alloc(strlen(prefix)+strlen(mesg)+1,sizeof(char));
    Rmesg = static_cast<char*>(Rheap);
    strcpy(Rmesg, prefix);
    strcat(Rmesg, mesg);
    return Rmesg;
}

