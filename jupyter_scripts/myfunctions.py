import numpy as np 
from netCDF4 import Dataset

ro2_names = ["c2h5o2","ch3o2", "hoch2o2", "hoch2ch2o2", "ch3choho2", \
             "ch3c(o)oo","hcoch2o2","hcoco3","hoch2co3"]

radical_names = ["ch3", "ch3o", "hco", "c2h5", "hoch2ch2o", \
                 "ch3co", "hcoco", "hochcho", "hoch2co"]
		

## Constants
kb = 1.38064852e-23

## Convert the string into a latex format chemical
## formula
def latex_name(string):
	latexname = ""

	for char in string:
		if char.isalpha():
			latexname = latexname + char.capitalize()
		if char.isdigit():
			latexname = latexname + "$_" + char + "$"

	return latexname

	
## Tracer molar masses 
mmol = {'ch4' : 16., 'ch3' : 15., 'ch3o2' : 47., \
		'ch3ooh' : 48., 'ch3oh' : 32., 'ch3o' : 31., \
		'hcho' : 30., 'hcooh' : 46., 'hoch2o2' : 63., \
		'hoch2oh' : 48., 'hoch2ooh' : 64., 'hco' : 29., \
		'c2h6' : 30., 'c2h5' : 29., 'c2h5o2' : 61. , \
		'c2h5ooh' : 62., 'c2h5oh' : 46., 'hoch2ch2o2' : 77., \
		'hoch2ch2o' : 61., 'ethgly' : 62., 'hyetho2h' : 78., \
		'ch3cho' : 44., 'ch2choh' : 44., 'ch3choho2' : 77., \
		'ch3cooh' : 60., 'ch3chohooh' : 78., 'ch3c(o)' : 43., \
		'ch3c(o)oo' : 75., 'ch3c(o)ooh' : 76., 'hcoch2o2' : 75., \
		'glyox' : 58., 'hcoco' : 57., 'hooch2cho' : 76., \
		'hoch2cho' : 60., 'hochcho' : 59., 'hoch2co' : 59., \
		'hoch2co3' : 91., 'hoch2co2h' : 76., 'hcoco2h' : 74., \
		'hcoco3h' : 90., 'hcoco3' : 89., 'hoch2co3h' : 92., \
		'co2' : 44., 'co' : 28., 'o' : 16., 'o1d' : 16., \
		'o2' : 32., 'o3' : 48., 'h' : 1., 'h2' : 2., 'oh' : 17., \
		'ho2' : 33., 'h2o2' : 34., 'n2' : 28., 'ar' : 40., \
		'h2o_vap' : 18., 'h2o_ice' : 18.,\
        'cl' : 35., 'cl2' : 37., 'hcl' : 36., 'hocl' : 52., \
       'clo' : 51., 'cloo' : 67., 'oclo' : 67., 'cl2o2' : 103.,\
        'ch3ocl' : 63., 'clco' : 63., 'clo3' : 83.5, 'hclo4' : 100.45,\
        'clo4' : 99.45\
        }
		
mmol_nonorg = {'co2' : 44., 'co' : 28., 'o' : 16., 'o1d' : 16., \
		'o2' : 32., 'o3' : 48., 'h' : 1., 'h2' : 2., 'oh' : 17., \
		'ho2' : 33., 'h2o2' : 34., \
		'h2o_vap' : 18., 'h2o_ice' : 18. }
		
mmol_org = {'ch4' : 16., 'ch3' : 15., 'ch3o2' : 47., \
		'ch3ooh' : 48., 'ch3oh' : 32., 'ch3o' : 31., \
		'hcho' : 30., 'hcooh' : 46., 'hoch2o2' : 63., \
		'hoch2oh' : 48., 'hoch2ooh' : 64., 'hco' : 29., \
		'c2h6' : 30., 'c2h5' : 29., 'c2h5o2' : 61. , \
		'c2h5ooh' : 62., 'c2h5oh' : 46., 'hoch2ch2o2' : 77., \
		'hoch2ch2o' : 61., 'ethgly' : 62., 'hyetho2h' : 78., \
		'ch3cho' : 44., 'ch2choh' : 44., 'ch3choho2' : 77., \
		'ch3cooh' : 60., 'ch3chohooh' : 78., 'ch3c(o)' : 43., \
		'ch3c(o)oo' : 75., 'ch3c(o)ooh' : 76., 'hcoch2o2' : 75., \
		'glyox' : 58., 'hcoco' : 57., 'hooch2cho' : 76., \
		'hoch2cho' : 60., 'hochcho' : 59., 'hoch2co' : 59., \
		'hoch2co3' : 91., 'hoch2co2h' : 76., 'hcoco2h' : 74., \
		'hcoco3h' : 90., 'hcoco3' : 89., 'hoch2co3h' : 92.}

## Atmospheric number density calculation 
def atmospheric_numdens(P,T):
	nd = P/(kb*T*1.e6)
	return nd 

## MMR -> VMR conversion 
def vmr(gas,mmr,mmean):
	gas = gas.replace("(","")
	gas = gas.replace(")","")
	return mmr*mmean/mmol[gas]
	
## MMR -> Number denisty conversion
def tracer_numdens(gas,nd,mmr,mmean):
	gas = gas.replace("(","")
	gas = gas.replace(")","")
	return vmr(gas,mmr,mmean)*nd
	
## Scan the rate coefficient titles and
## extract the production and loss rates
## of a specified gas
def rate_terms(ncdf,tracer):
    
    P_terms = {}
    
    L_terms = []
    
    for key in ncdf.variables.keys():
        var = ncdf[key]
        
        
        factor_1 = " " + tracer + " "
        factor_x = "*" + tracer + " "
        
        
        if ('title' in var.ncattrs() ):
            
            ## If the variable has save notation of rate coefficients...
            if ( "->" not in var.title):
                    continue
            else:       
                ## Find index of the arrow in the equation to
                ## split reactants from products
                arrow_index = var.title.find("->")
                products = var.title[arrow_index+2:]
                reactants = var.title[:arrow_index]
                
                #### Production
                ## If 1 mole of product comes from 1 mole of reactants...
                if ( factor_1 in products ):
                    P_terms[key] = 1.
                    continue
                ## If X mole of product comes from 1 mole of reactants...
                if ( factor_x in products):
                    
                    for p in products.split():
                        
                        if ( tracer in p ):
                            index_aster = p.find("*")
                            P_terms[key] = float(p[:index_aster])
                            continue
                #### Loss 
                for r in reactants.split():
                    if (r == tracer):
                        L_terms.append(key)
                            
    return P_terms, L_terms
							
					
	
def latex_title(string):
    
    if ( string == "o1d" ):
        title_latex = "O($^1$D)"
        return title_latex
    
    if ( string == "h2o_vap" ):
        title_latex = "H$_2$O Vapour"
        return title_latex 
    
    title_latex = ""
    
    if (string == "hoch2ooh + o2 + hv -> hcooh + ho2 + oh"):
        title_latex = "HOCH$_2$OOH + h$\\nu$ (+ O$_2$) $\longrightarrow$ HCOOH + HO$_2$ + OH"
        return title_latex
    
    for chunk in string.split():
     
        #### If the chunk is an operator
        #### we simply add to the latex
        #### string
        if ("->" == chunk):
            title_latex = title_latex + "$\longrightarrow$ "
            continue
        if ( "+" == chunk):
            title_latex = title_latex + chunk + " "
            continue
        
        if ("hv" == chunk):
            title_latex = title_latex + "h$\\nu$" + " "
            continue
            
        if ("o(1d)" == chunk):
            title_latex = title_latex + "O($^1$D) "
            continue
            
        if ("3ch2" == chunk):
            title_latex = title_latex + "$^3$CH$_2$ "
            continue
            
        if ("1ch2" == chunk):
            title_latex = title_latex + "$^1$CH$_2$ "
            continue
        
        #### Check if there's a factor
        if ( "*" in chunk):
            i = chunk.find("*")
        else:
            i = 0
        
        #### Seperate the factor from the gas
        gas = chunk[i:]
        
        latex_gas = ""
        for char in gas:
            
            if (char.isalpha()):
                latex_gas = latex_gas + char.upper()
            elif (char.isdigit()):
                latex_gas = latex_gas + "$_" + str(char) + "$"
            else:
                latex_gas = latex_gas + char

        if (i>0):
            latex_gas = chunk[:i] + latex_gas 
        
        title_latex = title_latex + latex_gas + " "
        
    return title_latex				
	
	
	
	