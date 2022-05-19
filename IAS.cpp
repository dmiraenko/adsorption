#include <iostream>
#include "functions.hpp"

int main(int argc, char* argv[]){

    pressure total_pressure = static_cast<pressure>(7 * pascal);

    dimless y[4] = {0.37, 0.13, 0.01, 0.49};
    conc_molpkg cs[4] = {32.47 * molpkg, 82.34 * molpkg, 584.43 * molpkg, 30.15 * molpkg};
    inv_pressure b[4] = {0.255 / (mega * pascal), 2.7682 / (mega * pascal), 97.7962 / (mega * pascal), 23.703 / (mega * pascal)};
    dimless t[4] = {0.777, 0.323, 0.134, 0.600};

    Unilan<conc_molpkg> f = Unilan<conc_molpkg> (cs[0], b[0], t[0]);
    pressure p = (pressure)(1e1 * atm);

    conc_molpkg z = 10 * molpkg;

    cout << f.surface_concetration(p) << " " 
        << f.henry_law_concentration(p) 
        << " " << f.reduced_spreading_pressure(p) 
        << " " << f.pure_component_pressure(z)
        << endl;



    // cout << name_format << engineering_prefix << static_cast<double>(langmuir(cs[0], p, b[0]) / cs[0]) << endl;
}