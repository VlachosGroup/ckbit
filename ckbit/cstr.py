"""
ckbit.cstr

"""

import pystan
from datetime import datetime
import numpy as np
import pickle
from hashlib import md5
import matplotlib.pyplot as plt
import pandas as pd
import arviz
from vunits.constants import R as R_conv 
from tabulate import tabulate

#Write the Stan code block for the CSTR model
def write_cstr_stan_code(runs, rxns, species, rxn_eqs, c_in, pH=False, \
                         R=0.0083144598, A0_up_lim=35, Ea_up_lim=350, \
                         priors=None):
    '''Writes Stan code that can be used for CSTR estimation
    
    Parameters
    ----------
        runs : int
            Number of experimental runs in dataset
        rxns : int
            Number of reactions in network
        species : int
            Number of species in dataset
        rxn_eqs : str
            Reaction equations
        c_in : list
            Inlet concentrations, each run is a list within the list
        pH : bool, optional
            Indicates whether pH is in the reaction equations and dataset, 
            Default is false
        R : float, optional
            Ideal gas constant for user selected units, Default is 0.0083144598
            which corresponds to kJ/mol/K, altering this value will necessitate
            altering Ea_up_lim and the priors for Ea accordingly
        A0_up_lim : float, optional
            Upper limit of A0 terms to be estimated, Default is 35
        Ea_up_lim : str, optional
            Upper limit of Ea terms to be estimated, Default is 350
        priors : list of str, optional
            User defined prior distributions, Must have appropriate format (see
            examples) in accordance with Stan, Default is None
            
    Returns
    -------
        code_cstr : str
            Code written in Stan syntax used for CSTR estimation
    '''

    #Building the data block
    data_block = '\n' + \
                'data {\n' + \
                '  int<lower=0> S;            // number of species measured\n' \
                + \
                '  int<lower=0> Rxns;         // number of reactions in network'
    for i in range(runs):
        n_runs = i+1
        n_times = '  int<lower=0> N{};           // number of measurement ' \
                   'times for run {}'.format(n_runs,n_runs)
        meas_times = '  vector<lower=0>[N{}] t{};    // measurement times ' \
                    'for run {}'.format(n_runs,n_runs,n_runs)
        temp = '  real<lower=0> T{};          // temperature of run {} ' \
               'in K'.format(n_runs,n_runs)
        concs = '  matrix<lower=0>[N{}, S] ci{}; // concentration ' \
                'measurements of species for run' \
                '{}'.format(n_runs,n_runs,n_runs)
        if pH:
            pH_val = '  real<lower=0> pH{};         // pH of run ' \
            '{}'.format(n_runs,n_runs)
            data_block = data_block + '\n' + n_times + '\n' + meas_times + \
                         '\n' + temp + '\n' + concs + '\n' + pH_val
        else:
            data_block = data_block + '\n' + n_times + '\n' + meas_times + \
                        '\n' + temp + '\n' + concs + '\n'
    data_block = data_block + '\n}\n'
    #Building the parameters block
    par_block = 'parameters {\n' + \
               '  vector<lower=0, upper={}>[Rxns] A0; // pre-exponential ' \
               'terms\n'.format(A0_up_lim) + \
               '  vector<lower=0, upper={}>[Rxns] Ea; // activation energy ' \
               'terms\n'.format(Ea_up_lim) + \
               '  real<lower=0> sigma;                // measurement error\n}\n'           
    #Building the transformed parameters block
    trans_par_block = 'transformed parameters{'
    for i in range(runs):
        n_runs = i+1
        for j in range(rxns):
            n_rxns = j+1
            k_term = '  real<lower=0> k_run{}_rxn{} = (10^A0[{}])*' \
                    'e()^(-Ea[{}]/({}*T{}))' \
                    ';'.format(n_runs,n_rxns,n_rxns,n_rxns,R,n_runs)
            trans_par_block = trans_par_block + '\n' + k_term
    trans_par_block = trans_par_block + '\n}\n'
    #Building the model block
    model_block = 'model {\n' + \
                 '  sigma ~ cauchy(0, 10);\n'
    for i in range(rxns):
        n_rxns = i + 1
        par_priors = '  A0[{}] ~ normal(10, 25);\n' \
                  '  Ea[{}] ~ normal(100, 250);\n'.format(n_rxns,n_rxns)
        model_block = model_block + par_priors
    if priors:
        for i in range(len(priors)):
            term = priors[i].split('~')[0]
            if model_block.find(term)==-1:
                    raise UserWarning('{}not found as a variable, cannot set' \
                                      ' its prior'.format(term)) 
            if priors[i].find('sigma')!=-1:
                model_block = model_block.replace('{}~ cauchy(0, 10)'.format(\
                                                  term),priors[i])
            if priors[i].find('A0')!=-1:
                model_block = model_block.replace('{}~ normal(10, 25)'.format(\
                                                  term),priors[i])
            if priors[i].find('Ea')!=-1:
                model_block = model_block.replace('{}~ normal(100, 250' \
                                                  ')'.format(term),priors[i])
    for i in range(runs):
        n_runs = i + 1
        model_block = model_block + '\n' + '  for (i in 1:N{}){{'.format(n_runs)
        for j in range(species):
            n_species = j + 1
            eq = rxn_eqs[j]
            eq = str.replace(eq, 't', 't{}[i]'.format(n_runs))
            eq = str.replace(eq, 'k', 'k_run{}_rxn'.format(n_runs))
            if pH:
                eq = str.replace(eq, 'pH', 'pH{}'.format(n_runs))
            else:
                eq = eq
            for m in range(species):
                n_species_eq = m + 1
                c_in_val = str(c_in[i][m])
                eq = str.replace(eq, 'c{}in'.format(n_species_eq), c_in_val)
                eq = str.replace(eq, 'c{}'.format(n_species_eq), \
                                      'ci{}[i,{}]'.format(n_runs,n_species_eq)) 
            
            model_block = model_block + '\n' + \
                         '    ci{}[i,{}] ~ normal('.format(n_runs, n_species) +\
                         eq + ', sigma);'
        model_block = model_block + '\n  }' 
    model_block = model_block + '\n}\n'   
    code_cstr = data_block+par_block+trans_par_block+model_block
    return code_cstr

#Experimental data import
def cstr_exp_data(filename, pH=False):
    '''Processes Excel file with CSTR data that can be used for CSTR estimation
    
    Parameters
    ----------
        filename : str
            Filename of Excel input file that contains the appropriately 
            formated CSTR data (see examples)
        pH : bool, optional
            Indicates whether pH is in the reaction equations and dataset, 
            Default is false
        
    Returns
    -------
        n_runs : int
            Number of experimental runs in dataset
        rxns : int
            Number of reactions in network
        S : int
            Number of species in dataset
        rxn_eqs : str
            Reaction equations
        c_in : list
            Inlet concentrations, each run is a list within the list
        cstr_data: dict
            Dictionary containing CSTR data inputs for Stan code
            
    '''
    file = pd.ExcelFile(filename)
    #Equation data processing
    equations = file.parse('Equations')
    S = len(equations.Species)
    rxn_eqs = equations.Equation
    rxn_counter = 0
    rxn_incomplete = True
    while rxn_incomplete:
        rxn_counter = rxn_counter + 1
        k_term = 'k{}'.format(rxn_counter)
        rxn_incomplete = np.any(equations.Equation.str.contains(k_term))
    rxns = rxn_counter-1
    cstr_data = {'S': S, 'Rxns': rxns}
    #Concentration data processing
    data = file.parse('Data')
    data = data.drop(data.index[0])
    n_runs = int(max(data.Run))
    c_in = []
    for i in range(n_runs):
        run = i + 1
        df = data[data['Run'] == run].reset_index(drop='True')
        cols = len(df.columns)
        c_in.append(np.array(df.iloc[0,(cols-S):(cols)].astype('float64')))
        df = df.drop(df.index[0])
        N = len(df.Run)
        t = df.Time
        T = df.Temp[1]+273
        c = np.array(df.iloc[:,(cols-S):(cols)].astype('float64'))
        newcstr_data = {'N{}'.format(run): N,
                       't{}'.format(run): t,
                       'T{}'.format(run): T,
                       'ci{}'.format(run): c,}
        if pH: 
            pH_val = df.pH[1]
            newcstr_data['pH{}'.format(run)] = pH_val
        cstr_data.update(newcstr_data)
    return n_runs, rxns, S, rxn_eqs, c_in, cstr_data
    

#Code to run CSTR MCMC Estimate
def MCMC(filename, model_name='cstr', pH=False, warmup=None, iters=2000, \
             chains=2, n_jobs=1, verbose=True, R_units='kJ/mol/K', \
             seed=None, \
             control={'adapt_delta':0.8, 'max_treedepth':20}, trace=True, \
             A0_up_lim=35, Ea_up_lim=350, priors=None, init_random=False, \
             A0_init=10, Ea_init=80, sigma_init=1):
    '''Bayesian inference using MCMC sampling for CSTR parameter estimation
    
    Parameters
    ----------
        filename : str
            Filename of Excel input file that contains the appropriately 
            formated CSTR data (see examples)
        model_name : str, optional
            Name of model, used for saving/loading compilied Stan code, Default
            is 'cstr'
        pH : bool, optional
            Indicates whether pH is in the reaction equations and dataset, 
            Default is false
        warmup : int, optional
            Number of warmup iterations for MCMC sampler for each chain, Must be
            less than the number of total iterations, Default is None, which 
            sets warmup equal to half of iters (the Stan default)
        iters : int, optional
            Number of total interations for MCMC sampler for each chain, Must be
            greater than the warmup, total number of samples useable for MCMC 
            inference will equal (chains*(iters-warmup)), Default is 2000
        chains : int, optional
            Number of chains for MCMC sampler, Default is 2
        n_jobs : int, optional
            Number of jobs to run in parallel for MCMC sampler, maximum is 
            number of cores the computer has, Default is 1
        verbose : bool, optional
            Flag to signal whether Stan intermediate output should be piped to
            terminal, Default is True
        R_units : str, optional
            Ideal gas constant units, Default is 'kJ/mol/K', note: if changing
            these units, the prior for Ea should be updated accordingly as
            well as Ea_up_lim
        seed : int, optional
            A positive integer used to seed the random number generation, use one seed even 
            when multiple chains are used since the other chain’s seeds are generated 
            from the first chain’s to avoid dependency among random number streams, set this
            seed for repeatable inference sampling runs, Default is np.random.randint(0, 1E9)
        control : dict, optional
            Dictionary of settings for MCMC sampler, Default is 
            {'adapt_delta':0.9999, 'max_treedepth':100}, more information at:
            https://mc-stan.org/docs/2_23/reference-manual/hmc-algorithm-parameters.html
        trace : bool, optional
            Flag to signal whether traceplots should be generated upon the run's
            completion, Default is True
        A0_up_lim : float, optional
            Upper limit of A0 terms to be estimated, Default is 35
        Ea_up_lim : str, optional
            Upper limit of Ea terms to be estimated, Default is 350
        priors : list of str, optional
            User defined prior distributions, Must have appropriate format (see
            examples) in accordance with Stan, Default is None
        init_random : bool, optional
            Flag to signal whether the initialization should be random or if it 
            should use user specified values, Default is False
        A0_init : float, optional
            Initialization point for the sampler for A0, Default is 10
        Ea_init : float, optional
            Initialization point for the sampler for Ea, Default is 80
        sigma_init : float, optional
            Initialization point for the sampler for sigma, Default is 1
            
    Returns
    -------
        fit : Stan object
            Object containing results of the MCMC run
        sample_vals : dict
            Dictionary of values collected by the MCMC sampler
    '''
    startTime = datetime.now()
    R = R_conv(R_units)
    #Process data
    n_runs, rxns, S, rxn_eqs, c_in, cstr_data = \
        cstr_exp_data(filename=filename, pH=pH)
    #Write stan code
    cstr_code = write_cstr_stan_code(runs=n_runs, rxns=rxns, species=S, \
                                     rxn_eqs=rxn_eqs, c_in=c_in, pH=pH, R=R, \
                                     A0_up_lim=A0_up_lim, Ea_up_lim=Ea_up_lim, \
                                     priors=priors)
    if warmup is None:
        warmup = int(iters/2)
    elif warmup>=iters:
        raise UserWarning('\nWarmup must be less than iters\nWarmup'\
                          'Entry:{}\nIters Entry:{}'.format(warmup, iters))
    #Compile stan model or open old one
    sm = StanModel_cache(model_code=cstr_code, model_name=model_name)
    #Write initialization list
    init_list = []
    for i in range(chains):
        A0 = []
        Ea = []
        for j in range(rxns):
            A0.append(A0_init)
            Ea.append(Ea_init)
        dict_init = {'A0':A0, 'Ea':Ea, 'sigma':sigma_init}
        init_list.append(dict_init)
    if init_random: init_list='random'
    #Run sampler
    if seed==None: seed=np.random.randint(0, 1E9)
    fit = sm.sampling(data=cstr_data, warmup=warmup, iter=iters, chains=chains,\
                      n_jobs=n_jobs, verbose=verbose, control=control, \
                      pars=['A0','Ea','sigma'], init=init_list, seed=seed)
    #Generate and print results
    print(fit)
    sample_vals = fit.extract(permuted=True)
    if trace: 
        traceplot = arviz.plot_trace(fit, compact=False)
        for i in range(rxns):
            traceplot[i,0].set_title('A0{}'.format(i+1))
            traceplot[i,1].set_title('A0{}'.format(i+1))
            traceplot[i+rxns,0].set_title('Ea{}'.format(i+1))
            traceplot[i+rxns,1].set_title('Ea{}'.format(i+1))
        plt.show()
    total_runtime = ((datetime.now() - startTime).total_seconds())/60
    print('Runtime (min): %.4f' % total_runtime)
    return fit, sample_vals

#Code to run CSTR VI Estimate
def VI(filename, model_name='cstr', pH=False, R_units='kJ/mol/K', priors=None,\
       A0_up_lim=35, Ea_up_lim=350, iters=2000000, algorithm='fullrank', \
       seed=None,\
       sample_file='./samples.csv', diagnostic_file='./diagnostics.csv',\
       verbose=True, grad_samples=1, elbo_samples=100, tol_rel_obj=0.01, \
       adapt_iter=50, eval_elbo=100, output_samples=10000, eta=0.2, \
       adapt_engaged=False, trace=True, init_random=False, \
       A0_init=10, Ea_init=80, sigma_init=1):
    '''Bayesian inference using VI for CSTR parameter estimation
    
    Parameters
    ----------
        filename : str
            Filename that contains the appropriately formated CSTR data (see 
            examples)
        model_name : str, optional
            Name of model, used for saving/loading compilied Stan code, Default
            is 'cstr'
        pH : bool, optional
            Indicates whether pH is in the reaction equations and dataset, 
            Default is false
        R_units : str, optional
            Ideal gas constant units, Default is 'kJ/mol/K', note: if changing
            these units, the prior for Ea should be updated accordingly as
            well as Ea_up_lim
        priors : list of str, optional
            User defined prior distributions, Must have appropriate format (see
            examples) in accordance with Stan, Default is None
        A0_up_lim : float, optional
            Upper limit of A0 terms to be estimated, Default is 35
        Ea_up_lim : str, optional
            Upper limit of Ea terms to be estimated, Default is 350
        iters : int, optional
            Maximum number of interations for VI sampler, Default is 2,000,000
        algorithm : str, optional
            Algorithm to use for VI, either 'meanfield' (for uncorrelated 
            posteriors) or 'fullrank' (for correlated posteriors), Default is
            'fullrank'
        seed : int, optional
            A positive integer used to seed the random number generation, 
            Default is np.random.randint(0, 1E9)            
        sample_file : str, optional
            Filename where the VI samples are saved to, Default is './samples.
            csv'
        diagnostic_file : str, optional
            Filename where the VI diagonstics are saved to, Default is './
            diagnostics.csv'
        verbose : bool, optional
            Flag to signal whether Stan intermediate output should be piped to
            terminal, Default is True
        grad_samples : int, optional
            Number of gradient evaluations to make to estimate gradient for VI
            solver, Default is 1
        elbo_samples : int, optional
            Number of elbo evaluations to make to estimate ELBO for VI solver, 
            Default is 100
        tol_rel_obj : float, optional
            Relative tolerance for VI solver, Default is 0.01
        adapt_iter : int, optional
            Number of iterations for adaptive tuning of eta, Default is 50
        eval_elbo : int, optional
            Number of iterations between ELBO evaluations for VI solver,
            Default is 100
        output_samples : int, optional
            Number of samples to draw from final approximation of posterior from
            VI solver, Default is 10,000
        eta : float, optional
            Positive, stepsize weighing parameter for VI solver, Ignored if 
            adapt_iter is True, Default is 0.2            
        adapt_engaged : 
            Flag to signal whether eta should be automatically tuned, Default is
            False
        trace : bool, optional
            Flag to signal whether traceplots should be generated upon the run's
            completion, Default is True
        init_random : bool, optional
            Flag to signal whether the initialization should be random or if it 
            should use user specified values, Default is False
        A0_init : float, optional
            Initialization point for the sampler for A0, Default is 10
        Ea_init : float, optional
            Initialization point for the sampler for Ea, Default is 80
        sigma_init : float, optional
            Initialization point for the sampler for sigma, Default is 1
 
    Returns
    -------
        fit : Stan object
            Stan object containing results of the VI run
        sample_vals : dict
            Dictionary of values collected by the VI sampler
    '''
    startTime = datetime.now()
    R = R_conv(R_units)
    #Process data
    n_runs, rxns, S, rxn_eqs, c_in, cstr_data = \
        cstr_exp_data(filename=filename, pH=pH)
    #Write stan code
    cstr_code = write_cstr_stan_code(runs=n_runs, rxns=rxns, species=S, \
                                     rxn_eqs=rxn_eqs, c_in=c_in, pH=pH, R=R, \
                                     A0_up_lim=A0_up_lim, Ea_up_lim=Ea_up_lim, \
                                     priors=priors)
    #Compile stan model or open old one
    sm = StanModel_cache(model_code=cstr_code, model_name=model_name)
    #Write initialization dictionary
    A0 = []
    Ea = []
    for j in range(rxns):
        A0.append(A0_init)
        Ea.append(Ea_init)
    dict_init = {'A0':A0, 'Ea':Ea, 'sigma':sigma_init}
    if init_random: dict_init='random'
    #Run variational inference
    if seed==None: seed=np.random.randint(0, 1E9)
    fit = sm.vb(data=cstr_data, algorithm=algorithm, iter=iters, \
                verbose=verbose, pars=['A0','Ea','sigma'], init=dict_init, \
                sample_file=sample_file, diagnostic_file=diagnostic_file, \
                grad_samples=grad_samples, elbo_samples=elbo_samples, \
                tol_rel_obj=tol_rel_obj, adapt_iter=adapt_iter, \
                adapt_engaged=adapt_engaged, eta=eta, \
                eval_elbo=eval_elbo, output_samples=output_samples, seed=seed)
    #Generate and print results
    sample_vals = fit['sampler_params']
    sample_names = fit['sampler_param_names']
    dict_vals = {}
    for i in range(2*rxns+1):
        dict_vals['{}'.format(sample_names[i])] = sample_vals[i]
    if trace: arviz.plot_trace(dict_vals)
    names = fit['sampler_param_names'][:-1]
    rows = len(names)
    data_table = []
    for i in range(0,rows):
        data_table.append([names[i], round(np.mean(sample_vals[i]),2), \
                         round(np.std(sample_vals[i]),2), \
                         round(np.quantile(sample_vals[i],0.025),2), \
                         round(np.quantile(sample_vals[i],0.25),2), \
                         round(np.quantile(sample_vals[i],0.5),2), \
                         round(np.quantile(sample_vals[i],0.75),2), \
                         round(np.quantile(sample_vals[i],0.975),2)])
    print(tabulate(data_table, headers=['', 'mean', 'sd', '2.5%', '25%',
                                 '50%', '75%', '97.5%']))
    with open(diagnostic_file, 'r') as f_ptr:
        lines = f_ptr.readlines()
    final_vals = lines[-1]
    iter_val, time_val, elbo_val = final_vals.split(',')
    if int(float(iter_val))==iters: print('The maximum number of iterations ' \
          'is reached! The algorithm may not have converged. Consider ' \
          'increasing the iters parameter by a factor of 10.')
    print('Check Convergence of ELBO plot to ensure ELBO converged corretly.' \
          ' The data points should approach and stabilize at a maximum'\
          'value, and there should be at least 10,000 iterations. If not ' \
          'converged, run again with a doubled eta value. Default eta value ' \
          'is 0.2 . It is recommended to run this twice with different ' \
          'random seed initializations and ensure the ' \
          'results are consistent.'.format(elbo_val))
    data = pd.read_csv(diagnostic_file ,skiprows=range(0,21), \
                       names=['iters','times','elbo'])
    iters67 = np.rint(0.67*len(data['elbo']))
    y_range = np.mean(data['elbo'][int(iters67):len(data['elbo'])])*2
    f, ax = plt.subplots(1)
    ax.scatter(data['iters'],data['elbo'])
    if y_range>0:
        ax.axes.set_ylim([0,y_range])
    elif y_range<0:
        ax.axes.set_ylim([y_range,0])
    ax.set_xlabel('Iterations')
    ax.set_ylabel('ELBO Value')
    ax.set_title('Convergence of ELBO')
    total_runtime = ((datetime.now() - startTime).total_seconds())/60
    print('Runtime (min): %.4f' % total_runtime)
    return fit, sample_vals

#Code to run CSTR MAP Estimate
def MAP(filename, model_name='cstr', pH=False, R_units='kJ/mol/K', priors=None,\
         A0_up_lim=35, Ea_up_lim=350, init_random=False, \
         seed=None, \
         verbose=True, A0_init=10, Ea_init=80, sigma_init=1):
    '''MAP estimation for CSTR parameter estimation
    
    Parameters
    ----------
        filename : str
            Filename of Excel input file that contains the appropriately 
            formated CSTR data (see examples)
        model_name : str, optional
            Name of model, used for saving/loading compilied Stan code, Default
            is 'cstr'
        pH : bool, optional
            Indicates whether pH is in the reaction equations and dataset, 
            Default is false
        R_units : str, optional
            Ideal gas constant units, Default is 'kJ/mol/K', note: if changing
            these units, the prior for Ea should be updated accordingly as
            well as Ea_up_lim
        priors : list of str, optional
            User defined prior distributions, Must have appropriate format (see
            examples), Default is None
        A0_up_lim : float, optional
            Upper limit of A0 terms to be estimated, Default is 35
        Ea_up_lim : str, optional
            Upper limit of Ea terms to be estimated, Default is 350
        init_random : bool, optional
            Flag to signal whether the initialization should be random or if it 
            should use user specified values, Default is False
        seed : int, optional
            A positive integer used to seed the random number generation, 
            Default is np.random.randint(0, 1E9)            
        verbose : bool, optional
            Flag to signal whether Stan intermediate output should be piped to
            terminal, Default is True
        A0_init : float, optional
            Initialization point for the sampler for A0, Default is 10
        Ea_init : float, optional
            Initialization point for the sampler for Ea, Default is 80
        sigma_init : float, optional
            Initialization point for the sampler for sigma, Default is 1
            
    Returns
    -------
        point_estimates : dict
            Dictionary containing values corresponding to modes of posterior
    '''
    startTime = datetime.now()
    R = R_conv(R_units)
    #Process data
    n_runs, rxns, S, rxn_eqs, c_in, cstr_data = \
        cstr_exp_data(filename=filename, pH=pH)
    #Write stan code
    cstr_code = write_cstr_stan_code(runs=n_runs, rxns=rxns, species=S, \
                                     rxn_eqs=rxn_eqs, c_in=c_in, pH=pH, R=R, \
                                     A0_up_lim=A0_up_lim, Ea_up_lim=Ea_up_lim, \
                                     priors=priors)
    #Compile stan model or open old one  
    sm = StanModel_cache(model_code=cstr_code, model_name=model_name)
    #Write initialization list
    init_list = []
    A0 = []
    Ea = []
    for j in range(rxns):
        A0.append(A0_init)
        Ea.append(Ea_init)
    dict_init = {'A0':A0, 'Ea':Ea, 'sigma':sigma_init}
    init_list.append(dict_init)
    if init_random: init_list='random'
    #Run point estimation
    if seed==None: seed=np.random.randint(0, 1E9)
    point_estimates = sm.optimizing(data=cstr_data, verbose=verbose,\
                                    init=init_list, seed=seed)
    #Generate and print results
    data_table = []
    for i in point_estimates:
        if i=='A0' or i=='Ea':
            if rxns==1:
                data_table.append([i,point_estimates[i]])
            else:
                count = 1
                for j in range(len(point_estimates[i])):
                    name = '{}[{}]'.format(i,count)
                    count = count + 1
                    val = round(float(point_estimates[i][j]),2)
                    data_table.append([name,val])
        elif i=='sigma':
            data_table.append([i,point_estimates[i]])
    print(tabulate(data_table, headers=['Parameter', 'Estimate']))     
    total_runtime = ((datetime.now() - startTime).total_seconds())/60
    print('Runtime (min): %.4f' % total_runtime)
    return point_estimates

#Saves/loads Stan models to avoid recompilation 
def StanModel_cache(model_code, model_name, **kwargs):
    '''Function for saving/loading compiled Stan code to avoid recompilation
    
    Parameters
    ----------
        model_code : str
            Stan code written in proper format
        model_name : str, optional
            Name of model, used for saving/loading compilied Stan code
        
    Returns
    -------
        sm : Stan model
            Stan object from pystan function StanModel
    '''
    code_hash = md5(model_code.encode('ascii')).hexdigest()
    if model_name is None:
        cache_fn = 'cached-model-{}.pkl'.format(code_hash)
    else:
        cache_fn = 'cached-{}-{}.pkl'.format(model_name, code_hash)
    try:
        f = open(cache_fn, 'rb')
    except:
        sm = pystan.StanModel(model_code=model_code)
        with open(cache_fn, 'wb') as f:
            pickle.dump(sm, f)
    else:
        sm = pickle.load(f)
        f.close()
        print("Using cached StanModel")
    return sm