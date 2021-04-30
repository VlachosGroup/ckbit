"""
ckbit.pfr

"""

import pystan
from datetime import datetime
import numpy as np
import pickle
from hashlib import md5
import pandas as pd
import arviz
from vunits.constants import R as R_conv 
import matplotlib.pyplot as plt
from tabulate import tabulate

#Write the Stan code block for the CSTR model
def write_pfr_stan_code(runs, rxns, species, rxn_eqs, R=0.0083144598,\
                     A0_up_lim=35, Ea_up_lim=350, pH=False, priors=None, \
                     rel_tol=5E-4, abs_tol=5E-4, max_num_steps=100):
    '''Writes Stan code that can be used for PFR estimation
    
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
        R : float, optional
            Ideal gas constant for user selected units, Default is 0.0083144598
            which corresponds to kJ/mol/K, altering this value will necessitate
            altering Ea_up_lim and the priors for Ea accordingly
        A0_up_lim : float, optional
            Upper limit of A0 terms to be estimated, Default is 35
        Ea_up_lim : str, optional
            Upper limit of Ea terms to be estimated, Default is 350
        pH : bool, optional
            Indicates whether pH is in the reaction equations and dataset, 
            Default is false
        priors : list of str, optional
            User defined prior distributions, Must have appropriate format (see
            examples) in accordance with Stan, Default is None
        rel_tol : float, optional
            Relative tolerance for Stan's ODE solver, Default is 5E-4
        abs_tol : float, optional
            Absolute tolerance for Stan's ODE solver, Default is 5E-4
        max_num_steps : int, optional
            Maximum number of steps for Stan's ODE solver, Default is 100
            
    Returns
    -------
        code_pfr : str
            Code written in Stan syntax used for PFR estimation
    '''
    #Building the functions block
    functions_block = '\n' + \
                     'functions {\n'
    for i in range(runs):
        n_runs = i+1
        formula = '  real[] dz_dt{}(real t{},      // time\n' \
                  '                real[] z{},    // concentration\n' \
                  '                real[] theta, // parameters\n' \
                  '                real[] x_r,   // unused data\n' \
                  '                int[] x_i) {{\n'.format(n_runs,n_runs,n_runs)
        functions_block = functions_block + formula
        for j in range(species):
            n_species = j+1
            species_term = '    real ci{}_{} = z{}[{}];\n'.format(n_runs,\
                                       n_species,n_runs,n_species)
            functions_block = functions_block + species_term
        functions_block = functions_block + '    real T{} = x_r[1];\n'.format(\
                                                       n_runs)
        if pH: 
            functions_block = functions_block + '    real pH{} = x_r[2];'\
                                                '\n'.format(n_runs)
        for k in range(rxns):
            n_rxns = k+1
            A0 = k+1
            Ea = (k+1)+rxns
            k_term = '    real k_run{}_rxn{} = (10^theta[{}])'.format(n_runs\
                                    ,n_rxns,A0) + \
                    '*e()^(-theta[{}]/({}*T{}));\n'.format(Ea,R,n_runs)
            functions_block = functions_block + k_term
        fun_returns = '    return {' 
        for l in range(species):
            n_species = l+1
            eq = rxn_eqs[l]
            eq = str.replace(eq, 'k', 'k_run{}_rxn'.format(n_runs))
            if pH:
                eq = str.replace(eq, 'pH', 'pH{}'.format(n_runs))
            else:
                eq = eq
            eq = str.replace(eq, 'c', 'ci{}_'.format(n_runs))
            eq = '    real dci{}_{}_dt = '.format(n_runs,n_species) + eq + ';\n'
            functions_block = functions_block + eq
            fun_returns = fun_returns + ' dci{}_{}_dt,'.format(n_runs,n_species)
        fun_returns = fun_returns[:-1] + ' };\n'
        functions_block = functions_block + fun_returns + '  }\n'
    functions_block = functions_block + '}\n'
    #Building the data block
    data_block = 'data {\n' + \
                '  int<lower=0> S;            // number of species measured\n'\
                + \
                '  int<lower=0> Rxns;         // number of reactions in network'
    for i in range(runs):
        n_runs = i+1
        n_times = '  int<lower=0> N{};           // number of measurement ' \
                   'times for run {}'.format(n_runs,n_runs)
        meas_times = '  real<lower=0> t{}[N{}];      // measurement times ' \
                    'for run {}'.format(n_runs,n_runs,n_runs)
        temp = '  real<lower=0> T{}[1];       // temperature of run {} ' \
               'in K'.format(n_runs,n_runs) 
        concs_init = '  real<lower=0> c0_{}[S];     // initial concentration ' \
                'measurements of species for run {}'.format(n_runs,n_runs)
        concs = '  real<lower=0> ci_{}[N{}, S]; // concentration ' \
                'measurements of species for run {}'.format(n_runs,n_runs,\
                                                 n_runs)
        if pH:
            pH_val = '  real<lower=0> pH{}[1];      // pH of run {}'.format(\
                                        n_runs,n_runs)
            data_block = data_block + '\n' + n_times + '\n' + meas_times + \
                        '\n' + temp + '\n' + concs_init + '\n' + concs + '\n' \
                        + pH_val
        else:
            data_block = data_block + '\n' + n_times + '\n' + meas_times + \
                        '\n' + temp + '\n' + concs_init + '\n' + concs + '\n'
    data_block = data_block + '\n}\n'
    #Building the parameters block
    par_block = 'parameters {\n' + \
               '  real<lower=0, upper={}> A0[Rxns];   // list of A0 ' \
               'parameters\n'.format(A0_up_lim) + \
               '  real<lower=0, upper={}> Ea[Rxns];   // list of Ea ' \
               'parameters\n'.format(Ea_up_lim) + \
               '  real<lower=0> sigma;       // measurement error\n}\n'
    #Building the transformed parameters block
    trans_par_block = 'transformed parameters {\n' + \
                    '  real theta[2*Rxns] = append_array(A0, Ea);\n' 
    for i in range(runs):
        n_runs = i+1
        if pH:
            ode_fun = '  real z{}[N{}, S]\n' \
                     '  = integrate_ode_bdf(dz_dt{}, c0_{}, 0, t{}, theta,\n' \
                     '                     append_array(T{}, pH{}),' \
                     ' rep_array(0, 0),\n' \
                     '                     {}, {}, {});' \
                     '\n'.format(n_runs,n_runs,n_runs,n_runs,n_runs,n_runs,\
                                 n_runs,rel_tol,abs_tol,max_num_steps)
        else:
            ode_fun = '  real z{}[N{}, S]\n' \
                     '  = integrate_ode_bdf(dz_dt{}, c0_{}, 0, t{}, theta,\n' \
                     '                     T{}, rep_array(0, 0),\n' \
                     '                     {},{},{});' \
                     '\n'.format(n_runs,n_runs,n_runs,n_runs,n_runs,n_runs,\
                                 rel_tol,abs_tol,max_num_steps)
        trans_par_block = trans_par_block + ode_fun        
    trans_par_block = trans_par_block + '}\n'
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
        model_block = model_block + '  for (i in 1:N{}){{\n'.format(n_runs)
        model_block = model_block + '    ci_{}[i, ] ~ normal(z{}[' \
                                  'i, ], sigma);\n'.format(n_runs,n_runs)
        
        model_block = model_block + '  }\n' 
    model_block = model_block + '}\n'   
    code_pfr = functions_block+data_block+par_block+trans_par_block+model_block
    return code_pfr

#Experimental data import
def pfr_exp_data(filename, pH=False):
    '''Processes Excel file with PFR data that can be used for PFR estimation
    
    Parameters
    ----------
        filename : str
            Filename of Excel input file that contains the appropriately 
            formated PFR data (see examples)
        pH : bool, optional
            Indicates whether pH is in the reaction equations and dataset, 
            Default is false
        
    Returns
    -------
        n_runs : int
            Number of experimental runs in dataset
        rxns : int
            Number of reactions in network
        sp : int
            Number of species in dataset
        rxn_eqs : str
            Reaction equations
        pfr_data: dict
            Dictionary containing PFR data inputs for Stan code
            
    '''
    file = pd.ExcelFile(filename)
    #Equation data processing
    equations = file.parse('Equations')
    sp = len(equations.Species)
    rxn_eqs = equations.Equation
    rxn_counter = 0
    rxn_incomplete = True
    while rxn_incomplete:
        rxn_counter = rxn_counter + 1
        k_term = 'k{}'.format(rxn_counter)
        rxn_incomplete = np.any(equations.Equation.str.contains(k_term))
    rxns = rxn_counter-1
    pfr_data = {'S': sp, 'Rxns': rxns}
    #Concentration data processing
    data = file.parse('Data')
    data = data.drop(data.index[0])
    n_runs = int(max(data.Run))
    for i in range(n_runs):
        run = i + 1
        df = data[data['Run'] == run].reset_index(drop='True')
        cols = len(df.columns)
        cIn = (np.array(df.iloc[0,(cols-sp):(cols)].astype('float64')))
        df = df.drop(df.index[0])
        N = len(df.Run)
        t = df.Time
        T = df.Temp[1]+273
        c = np.array(df.iloc[:,(cols-sp):(cols)].astype('float64'))
        new_pfr_data = {'N{}'.format(run): N,
                      't{}'.format(run): t,
                      'T{}'.format(run): [T],
                      'c0_{}'.format(run): cIn,
                      'ci_{}'.format(run): c,}
        if pH: 
            pH_val = df.pH[1]
            new_pfr_data['pH{}'.format(run)] = [pH_val]
        pfr_data.update(new_pfr_data)
    return n_runs, rxns, sp, rxn_eqs, pfr_data

#Code to run PFR MCMC Estimate
def MCMC(filename, model_name='pfr', pH=False, R_units='kJ/mol/K', priors=None,\
         A0_up_lim=35, Ea_up_lim=350, warmup=None, iters=2000, chains=2, \
         n_jobs=1, verbose=True, trace=True, init_random=False,\
         seed=None,\
         control={'adapt_delta':0.8, 'max_treedepth':20}, rel_tol=5E-4, \
         abs_tol=5E-4, max_num_steps=100, A0_init=10, Ea_init=80, sigma_init=1):
    '''Bayesian inference using MCMC sampling for PFR parameter estimation
    
    Parameters
    ----------
        filename : str
            Filename of Excel input file that contains the appropriately 
            formated PFR data (see examples)
        model_name : str, optional
            Name of model, used for saving/loading compilied Stan code, Default
            is 'pfr'
        pH : bool, optional
            Flag indicating whether pH is in the reaction equations and dataset, 
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
        trace : bool, optional
            Flag to signal whether traceplots should be generated upon the run's
            completion, Default is True
        init_random : bool, optional
            Flag to signal whether the initialization should be random or if it 
            should use user specified values, Default is False
        seed : int, optional
            A positive integer used to seed the random number generation, use one seed even 
            when multiple chains are used since the other chain’s seeds are generated 
            from the first chain’s to avoid dependency among random number streams, set this
            seed for repeatable inference sampling runs, Default is np.random.randint(0, 1E9)
        control : dict, optional
            Dictionary of settings for MCMC sampler, Default is 
            {'adapt_delta':0.9999, 'max_treedepth':100}, more information at:
            https://mc-stan.org/docs/2_23/reference-manual/hmc-algorithm-parameters.html
        int_tol : float, optional
            Relative tolerance for Stan's ODE solver, Default is 5E-4
        abs_tol : float, optional
            Absolute tolerance for Stan's ODE solver, Default is 5E-4
        max_num_steps : int, optional
            Maximum number of steps for Stan's ODE solver, Default is 100
        A0_init : float, optional
            Initialization point for the sampler for A0, Default is 10
        Ea_init : float, optional
            Initialization point for the sampler for Ea, Default is 80
        sigma_init : float, optional
            Initialization point for the sampler for sigma, Default is 1
            
    Returns
    -------
        fit : Stan object
            Stan object containing results of the MCMC run
        sample_vals : dict
            Dictionary of values collected by the MCMC sampler
    '''
    startTime = datetime.now()
    R = R_conv(R_units)
    #Process data
    n_runs, rxns, sp, rxn_eqs, pfr_data = pfr_exp_data(filename=filename, pH=pH)
    #Write stan code
    pfr_code = write_pfr_stan_code(R=R, runs=n_runs, rxns=rxns, species=sp, \
                               rxn_eqs=rxn_eqs, pH=pH, priors=priors, \
                               A0_up_lim=A0_up_lim, Ea_up_lim=Ea_up_lim, \
                               rel_tol=rel_tol, abs_tol=abs_tol, \
                               max_num_steps=max_num_steps)
    if warmup is None:
        warmup = int(iters/2)
    elif warmup>=iters:
        raise UserWarning('\nWarmup must be less than iters\nWarmup'\
                          'Entry:{}\nIters Entry:{}'.format(warmup, iters))
    #Compile stan model or open old one
    sm = StanModel_cache(model_code=pfr_code, model_name=model_name)
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
    fit = sm.sampling(data=pfr_data, warmup=warmup, iter=iters, chains=chains, \
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

#Code to run PFR VI Estimate
def VI(filename, model_name='pfr', pH=False, R_units='kJ/mol/K', priors=None,\
       A0_up_lim=35, Ea_up_lim=350, iters=2000000, algorithm='fullrank', \
       verbose=True, rel_tol=5E-4, abs_tol=5E-4, max_num_steps=100, \
       seed=None, \
       sample_file='./samples.csv', diagnostic_file='./diagnostics.csv',\
       grad_samples=1, elbo_samples=100, tol_rel_obj=0.01, adapt_iter=50, \
       eval_elbo=100, output_samples=10000, eta=0.2, \
       adapt_engaged=False, trace=True, init_random=False, \
       A0_init=10, Ea_init=80, sigma_init=1):
    '''Bayesian inference using VI for PFR parameter estimation
    
    Parameters
    ----------
        filename : str
            Filename of Excel input file that contains the appropriately 
            formated PFR data (see examples)
        model_name : str, optional
            Name of model, used for saving/loading compilied Stan code, Default
            is 'pfr'
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
        verbose : bool, optional
            Flag to signal whether Stan intermediate output should be piped to
            terminal, Default is True
        rel_tol : float, optional
            Relative tolerance for Stan's ODE solver, Default is 5E-4
        abs_tol : float, optional
            Absolute tolerance for Stan's ODE solver, Default is 5E-4
        max_num_steps : int, optional
            Maximum number of steps for Stan's ODE solver, Default is 100
        seed : int, optional
            A positive integer used to seed the random number generation, 
            Default is np.random.randint(0, 1E9)            
        sample_file : str, optional
            Filename where the VI samples are saved to, Default is './samples.
            csv'
        diagnostic_file : str, optional
            Filename where the VI diagonstics are saved to, Default is './
            diagnostics.csv'
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
    #Process Data
    n_runs, rxns, sp, rxn_eqs, pfr_data = pfr_exp_data(filename=filename, pH=pH)
    #Write stan code
    pfr_code = write_pfr_stan_code(R=R, runs=n_runs, rxns=rxns, species=sp, \
                               rxn_eqs=rxn_eqs, pH=pH, priors=priors, \
                               A0_up_lim=A0_up_lim, Ea_up_lim=Ea_up_lim, \
                               rel_tol=rel_tol, abs_tol=abs_tol, \
                               max_num_steps=max_num_steps)
    #Compile stan model or open old one
    sm = StanModel_cache(model_code=pfr_code, model_name=model_name)
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
    fit = sm.vb(data=pfr_data, algorithm=algorithm, iter=iters, \
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

#Code to run PFR MAP Estimate
def MAP(filename, model_name='pfr', pH=False, R_units='kJ/mol/K', priors=None,\
         A0_up_lim=35, Ea_up_lim=350, init_random=False, \
         seed=None, \
         verbose=True, rel_tol=5E-4, abs_tol=5E-4, max_num_steps=100, \
         A0_init=10, Ea_init=80, sigma_init=1):
    '''MAP estimation for PFR parameter estimation
    
    Parameters
    ----------
        filename : str
            Filename of Excel input file that contains the appropriately 
            formated PFR data (see examples)
        model_name : str, optional
            Name of model, used for saving/loading compilied Stan code, Default
            is 'pfr'
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
        rel_tol : float, optional
            Relative tolerance for Stan's ODE solver, Default is 5E-4
        abs_tol : float, optional
            Absolute tolerance for Stan's ODE solver, Default is 5E-4
        max_num_steps : int, optional
            Maximum number of steps for Stan's ODE solver, Default is 100
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
    n_runs, rxns, sp, rxn_eqs, pfr_data = pfr_exp_data(filename=filename, pH=pH)
    #Generate stan code
    pfr_code = write_pfr_stan_code(R=R, runs=n_runs, rxns=rxns, species=sp, \
                               rxn_eqs=rxn_eqs, pH=pH, priors=priors, \
                               A0_up_lim=A0_up_lim, Ea_up_lim=Ea_up_lim, \
                               rel_tol=rel_tol, abs_tol=abs_tol, \
                               max_num_steps=max_num_steps)
    #Compile stan model or open old one
    sm = StanModel_cache(model_code=pfr_code, model_name=model_name)
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
    point_estimates = sm.optimizing(data=pfr_data, verbose=verbose,\
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