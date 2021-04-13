#ifndef FORALLOPERATIONS_H
#define FORALLOPERATIONS_H
#pragma once //just include once
#include<iostream>
#include<vector>
#include<cassert>
using std::vector;

template<typename T>
class forAllOperations
{
public:
    forAllOperations<T>();
    virtual ~forAllOperations<T>();
    
    vector<T> temp1dvector;
    vector<vector<T> > temp2dvector;
    vector<vector<vector<T> > > temp3dvector;

};
// defining macros

#define forAll(temp3dvector) \
for(unsigned int i = 0; i < temp3dvector.size(); i++) \
	for(unsigned int j = 0; j< temp3dvector[i].size(); j++) \
		for(unsigned int k = 0; k< temp3dvector[i][j].size(); k++)

#define forAllInternal(temp3dvector) \
for(unsigned int i = 1; i < temp3dvector.size()-1; i++) \
	for(unsigned int j = 1; j< temp3dvector[i].size()-1; j++) \
		for(unsigned int k = 1; k< temp3dvector[i][j].size() - 1; k++)

#define forAllInternalUCVs(temp3dvector) \
for(unsigned int i = 1; i<temp3dvector.size() -2; i++) \
	for(unsigned int j = 1; j<temp3dvector[i].size() -1; j++) \
		for(unsigned int k = 1; k<temp3dvector[i][j].size() - 1; k++)

#define forAllInternalVCVs(temp3dvector) \
for(unsigned int i = 1; i<temp3dvector.size() - 1; i++) \
	for(unsigned int j = 1; j<temp3dvector[i].size() -2; j++) \
		for(unsigned int k = 1; k<temp3dvector[i][j].size() - 1; k++)

#define forAllInternalWCVs(temp3dvector) \
for(unsigned int i = 1; i<temp3dvector.size() - 1; i++) \
	for(unsigned int j = 1; j<temp3dvector[i].size() - 1; j++) \
		for(unsigned int k = 1; k<temp3dvector[i][j].size() - 2; k++)


//for inlet boundary
#define forBottomBoundary(temp3dvector) \
for(unsigned int i = 1; i<temp3dvector.size()-1; i++) \
	for(unsigned int j = 1; j<temp3dvector[i].size() - 1; j++) \
		for(unsigned int  k = 0; k < 1; k++)

//for outlet boundary
#define forTopBoundary(temp3dvector) \
for(unsigned int i = 1; i<temp3dvector.size()-1; i++) \
	for(unsigned int j = 1; j<temp3dvector[i].size() - 1; j++) \
		for(unsigned int  k = temp3dvector[i][j].size() - 1; k < temp3dvector[i][j].size(); k++)


//wall
#define forNorthBoundary(temp3dvector) \
for(unsigned int i = 1; i<temp3dvector.size()-1; i++) \
	for(unsigned int j = temp3dvector[i].size() - 1; j<temp3dvector[i].size(); j++) \
		for(unsigned int  k = 1; k < temp3dvector[i][j].size() - 1; k++)

//Wall
#define forSouthBoundary(temp3dvector) \
for(unsigned int i = 1; i<temp3dvector.size()-1; i++) \
	for(unsigned int j = 0; j<1; j++) \
		for(unsigned int  k = 1; k < temp3dvector[i][j].size() - 1; k++)


//Wall
#define forEastBoundary(temp3dvector) \
for(unsigned int i = temp3dvector.size() - 1; i<temp3dvector.size(); i++) \
	for(unsigned int j = 1; j<temp3dvector[i].size() - 1; j++) \
		for(unsigned int  k = 1; k < temp3dvector[i][j].size() - 1; k++)


//Wall
#define forWestBoundary(temp3dvector) \
for(unsigned int i = 0; i<1; i++) \
	for(unsigned int j = 1; j<temp3dvector[i].size() - 1; j++) \
		for(unsigned int  k = 1; k < temp3dvector[i][j].size() - 1; k++)


#endif
