"""
ckbit.app_ea

"""

import pystan
from datetime import datetime
import numpy as np
import pickle
from hashlib import md5
import pandas as pd
import matplotlib.pyplot as plt
import arviz
from vunits.constants import R as R_conv 
from tabulate import tabulate

def write_app_ea_stan_code(Ea_up_lim=350, priors=None):
    '''Writes Stan code used for apparent activation energy estimation
    
    Parameters
    ----------
        Ea_up_lim : str, optional
            Upper limit of Ea terms to be estimated, Default is 350
        priors : list of str, optional
            User defined prior distributions, Must have appropriate format (see
            examples) in accordance with Stan, Default is None
       
    Returns
    -------
        code_app_ea : str
            Code written in Stan syntax used for apparent activation
            energy estimation
    '''
    #Building the data block
    data_block = 'data {\n' \
                 '  int<lower=0> N;\n' \
                 '  vector<upper=0>[N] x;\n' \
                 '  vector[N] y;\n' \
                 '}\n'
    #Building the parameters block
    par_block = 'parameters {\n' \
                '  real intercept;          // intercept of best ' \
                'fit line\n' + \
                '  real<lower=0, upper={}>'.format(Ea_up_lim) + \
                ' app_ea;   // slope of best fit line\n' \
               '  real<lower=0> sigma;               // measurement error\n' \
               '}\n'
    #Building the model block
    model_block = 'model {\n' \
                 '  sigma ~ cauchy(0, 10);\n' \
                 '  intercept ~ normal(10,100);\n' \
                 '  app_ea ~ normal(100,250);\n' \
                 '  y ~ normal(intercept + app_ea * x, sigma);' \
                 '}\n'
    if priors:
        for i in range(len(priors)):
            term = priors[i].split('~')[0]
            if model_block.find(term)==-1:
                    raise UserWarning('{}not found as a variable, cannot set' \
                                      ' its prior'.format(term)) 
            if priors[i].find('sigma')!=-1:
                model_block = model_block.replace('{}~ cauchy(0, 10)'.format(\
                                                  term),priors[i])
            if priors[i].find('intercept')!=-1:
                model_block = model_block.replace('{}~ normal(10,25)'.format(\
                                                  term),priors[i])
            if priors[i].find('app_ea')!=-1:
                model_block = model_block.replace('{}~ normal(100,250)'.format(\
                                                  term),priors[i])
    code_app_ea = data_block+par_block+model_block
    return code_app_ea

def app_ea_exp_data(filename, R_units = 'kJ/mol/K'):
    '''Processes Excel file with apparent activation energy data
    
    Parameters
    ----------
        filename : str
            Filename of Excel input file that contains the appropriately 
            formated apparent activation energy data (see examples)
        R_units : str, optional
            Ideal gas constant units, Default is 'kJ/mol/K', note: if changing
            these units, the prior for app_ea should be updated accordingly
        
    Returns
    -------
        app_ea_data: dict
            Dictionary containing apparent activation energy data inputs for 
            inputs into Stan code
            
    '''    
    #Experimental data import
    file = pd.ExcelFile(filename)
    data = file.parse('Data')
    temps = data.Temperature
    rates = data.Rate
    R = R_conv(R_units)
    #Data processing
    if ~(data.Temperature>0).all():
        raise Exception('Input Data Error: Temperature values must be ' \
                        'greater than 0 K') 
    scaled_temps = (-1/R)*(1/np.array(temps))
    lnRates = np.array(np.log(np.absolute(rates)))
    app_ea_data = {'N': np.size(scaled_temps),
                 'x': scaled_temps,
                 'y': lnRates}
    return app_ea_data

#Code to run apparent activation energy MCMC estimate
def MCMC(filename, model_name='app_ea', R_units='kJ/mol/K', priors=None,\
         Ea_up_lim=350, warmup=None, iters=5000, chains=2, n_jobs=1, \
         verbose=True, trace=True, init_random=False,\
         seed=None,\
         control={'adapt_delta':0.9999, 'max_treedepth':100}, int_init=10, \
         Ea_init=80, sigma_init=1):
    '''Bayesian inference using MCMC sampling for apparent Ea estimation
    
    Parameters
    ----------
        filename : str
            Filename of Excel input file that contains the appropriately 
            formated Ea data (see examples)
        model_name : str, optional
            Name of model, used for saving/loading compilied Stan code, Default
            is 'app_ea'
        R_units : str, optional
            Ideal gas constant units, Default is 'kJ/mol/K', note: if changing
            these units, the prior for app_ea should be updated accordingly
        priors : list of str, optional
            User defined prior distributions, Must have appropriate format (see
            examples) in accordance with Stan, Default is None
        Ea_up_lim : str, optional
            Upper limit of Ea terms to be estimated, Default is 350, note: if
            R is changed, this will likely need to be changed too
        warmup : int, optional
            Number of warmup iterations for MCMC sampler for each chain, Must be
            less than the number of total iterations, Default is None, which 
            sets warmup equal to half of iters (the Stan default)
        iters : int, optional
            Number of total interations for MCMC sampler for each chain, Must be
            greater than the warmup, total number of samples useable for MCMC 
            inference will equal (chains*(iters-warmup)), Default is 5000
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
        int_init : float, optional
            Initialization point for the sampler for intercept, Default is 10
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
    #Process data
    app_ea_data = app_ea_exp_data(filename=filename, R_units=R_units)
    #Write stan code
    app_ea_code = write_app_ea_stan_code(Ea_up_lim=Ea_up_lim, priors=priors)
    if warmup is None:
        warmup = int(iters/2)
    elif warmup>=iters:
        raise UserWarning('\nWarmup must be less than iters\nWarmup'\
                          'Entry:{}\nIters Entry:{}'.format(warmup, iters))
    #Compile stan model or open old one
    sm = StanModel_cache(model_code=app_ea_code, model_name=model_name)
    #Write initialization list
    init_list = []
    for i in range(chains):
        dict_init = {'intercept':int_init, 'app_ea':Ea_init, 'sigma':sigma_init}
        init_list.append(dict_init)
    if init_random: init_list='random'
    #Run sampler
    if seed==None: seed=np.random.randint(0, 1E9)
    fit = sm.sampling(data=app_ea_data, warmup=warmup, iter=iters, \
                      chains=chains, n_jobs=n_jobs, verbose=verbose, \
                      control=control, pars=['intercept','app_ea','sigma'], \
                      init=init_list, seed=seed)
    #Generate and print results
    print(fit)
    sample_vals = fit.extract(permuted=True)
    if trace: arviz.plot_trace(fit)
    total_runtime = ((datetime.now() - startTime).total_seconds())/60
    print('Runtime (min): %.4f' % total_runtime)
    return fit, sample_vals

#Code to run apparent activation energy VI estimate
def VI(filename, model_name='app_ea', R_units='kJ/mol/K', priors=None,\
         Ea_up_lim=350, iters=2000000, algorithm='fullrank', \
         verbose=True, seed=None,\
         sample_file='./samples.csv', diagnostic_file='./diagnostics.csv',\
         grad_samples=1, elbo_samples=100, tol_rel_obj=0.01, adapt_iter=50, \
         eval_elbo=100, output_samples=10000, eta=0.2, \
         adapt_engaged=False, trace=True, init_random=False,\
         int_init=10, Ea_init=80, sigma_init=1):
    '''Bayesian inference using VI for apparent Ea estimation
    
    Parameters
    ----------
        filename : str
            Filename of Excel input file that contains the appropriately 
            formated Ea data (see examples)
        model_name : str, optional
            Name of model, used for saving/loading compilied Stan code, Default
            is 'app_ea'
        R_units : str, optional
            Ideal gas constant units, Default is 'kJ/mol/K', note: if changing
            these units, the prior for app_ea should be updated accordingly
        priors : list of str, optional
            User defined prior distributions, Must have appropriate format (see
            examples) in accordance with Stan, Default is None
        Ea_up_lim : str, optional
            Upper limit of Ea terms to be estimated, Default is 350, this may
            need to be changed if changing the units of R
        iters : int, optional
            Maximum number of iterations of the variational density to achieve
            minimizing the ELBO before the VI algorithm terminates, Default is
            2,000,000
        algorithm : str, optional
            Algorithm to use for VI, either 'meanfield' (for uncorrelated 
            posteriors) or 'fullrank' (for correlated posteriors), Default is
            'fullrank' 
        verbose : bool, optional
            Flag to signal whether Stan intermediate output should be piped to
            terminal, Default is True
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
            Number of ELBO evaluations to make to estimate ELBO for VI solver, 
            Default is 100
        tol_rel_obj : float, optional
            Relative tolerance convergence criteria for median and mean of the 
            change in the ELBO for VI solver, Default is 0.01
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
            adapt_engaged is True, Default is 0.2  
        adapt_engaged : 
            Flag to signal whether eta should be automatically tuned, Default is
            False
        trace : bool, optional
            Flag to signal whether traceplots should be generated upon the run's
            completion, Default is True
        int_init : float, optional
            Initialization point for the sampler for intercept, Default is 10
        Ea_init : float, optional
            Initialization point for the sampler for Ea, Default is 80
        sigma_init : float, optional
            Initialization point for the sampler for sigma, Default is 1
        init_random : bool, optional
            Flag to signal whether the initialization should be random or if it 
            should use user specified values, Default is False

    Returns
    -------
        fit : Stan object
            Stan object containing results of the VI run
        sample_vals : dict
            Dictionary of values collected by the VI sampler
    '''
    startTime = datetime.now()
    #Process data
    app_ea_data = app_ea_exp_data(filename=filename, R_units=R_units)
    #Write stan code
    app_ea_code = write_app_ea_stan_code(Ea_up_lim=Ea_up_lim, priors=priors)
    #Compile stan model or open old one
    sm = StanModel_cache(model_code=app_ea_code, model_name=model_name)
    #Write initialization list
    dict_init = {'intercept':int_init, 'app_ea':Ea_init, 'sigma':sigma_init}
    if init_random: dict_init='random'
    #Run VI estimation
    if seed==None: seed=np.random.randint(0, 1E9)
    fit = sm.vb(data=app_ea_data, algorithm=algorithm, iter=iters, \
                verbose=verbose, seed=seed, init=dict_init,\
                sample_file=sample_file, diagnostic_file=diagnostic_file, \
                grad_samples=grad_samples, elbo_samples=elbo_samples, \
                tol_rel_obj=tol_rel_obj, adapt_iter=adapt_iter, \
                adapt_engaged=adapt_engaged, eta=eta, \
                eval_elbo=eval_elbo, output_samples=output_samples)
    #Generate and print results
    sample_vals = fit['sampler_params']
    sample_names = fit['sampler_param_names']
    dict_vals = {}
    for i in range(len(sample_vals)):
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

#Code to run apparent activation energy MAP estimate
def MAP(filename, model_name='app_ea', R_units='kJ/mol/K', priors=None,\
         Ea_up_lim=350, verbose=True, init_random=False,\
         seed=None,\
         int_init=10, Ea_init=80, sigma_init=1):
    '''MAP estimation for apparent Ea estimation
    
    Parameters
    ----------
        filename : str
            Filename that contains the appropriately formated Ea data (see 
            examples)
        model_name : str, optional
            Name of model, used for saving/loading compilied Stan code, Default
            is 'app_ea'
        R_units : str, optional
            Ideal gas constant units, Default is 'kJ/mol/K', note: if changing
            these units, the prior for app_ea should be updated accordingly
        priors : list of str, optional
            User defined prior distributions, Must have appropriate format (see
            examples) in accordance with Stan, Default is None
        Ea_up_lim : str, optional
            Upper limit of Ea terms to be estimated, Default is 350
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
            A positive integer used to seed the random number generation
            Default is np.random.randint(0, 1E9)
        int_init : float, optional
            Initialization point for the sampler for intercept, Default is 10
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
    #Process data
    app_ea_data = app_ea_exp_data(filename=filename, R_units=R_units)
    #Write stan code
    app_ea_code = write_app_ea_stan_code(Ea_up_lim=Ea_up_lim, priors=priors)
   #Compile stan model or open old one
    sm = StanModel_cache(model_code=app_ea_code, model_name=model_name)
    #Write initialization list
    init_list = [{'intercept':int_init, 'app_ea':Ea_init, 'sigma':sigma_init}]
    if init_random: init_list='random'
    #Run point estimation
    if seed==None: seed=np.random.randint(0, 1E9)
    point_estimates = sm.optimizing(data=app_ea_data, verbose=verbose, \
                                    init=init_list, seed=seed)
    #Generate and print results
    data_table = []
    for i in point_estimates:
        data_table.append([i,round(float(point_estimates[i]),2)])
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
