classdef SB_Params
    %SB_PARAMS params for a set of spin-boson simulations
    %   this will contain info like the number of simulations, various
    %   parameters, and a list of independent parameters that are varied
    
    properties
        lx
        ly
        a
        tc
        cutoff_type
        tcj
        gamma
        beta
        h
        alpha
        delta
        v
        omega_c
        n_sims
    end
    
    methods
        function obj = SB_Params(infile, N)
            obj.n_sims = N;
            fid = fopen(infile);
            tline = fgetl(fid);%#Lattice Parameters
            tline = fgetl(fid);%#dimension
            tline = fgetl(fid);%#lattice type
            
            tline = fgetl(fid);%#lengths
            C = textscan(tline, '%f');
            obj.lx = C{1}(1);
            obj.ly = C{1}(2);
            
            tline = fgetl(fid);%#spacings
            C = textscan(tline, '%f');
            obj.a = C{1}(1);
            obj.tc = C{1}(2);
            
            tline = fgetl(fid);%#blank
            tline = fgetl(fid);%#Model Parameters
            tline = fgetl(fid);%#cutoff type
            C = textscan(tline, '%s');
            obj.cutoff_type = C{1}{1};
            
            tline = fgetl(fid);%#J_couples
            C = textscan(tline, '%f');
            obj.tcj = C{1}(1);
            obj.gamma = C{1}(2);
            
            tline = fgetl(fid);%#beta
            C = textscan(tline, '%f');
            obj.beta = C{1}(1);
            
            tline = fgetl(fid);%#h
            C = textscan(tline, '%f');
            obj.h = C{1}(1);
            
            tline = fgetl(fid);%#blank
            tline = fgetl(fid);%#Algorithm Parameters
            tline = fgetl(fid);%#algorithm
            tline = fgetl(fid);%#random seed
            tline = fgetl(fid);%#blank
            tline = fgetl(fid);%#Simulation Parameters
            tline = fgetl(fid);%#equilibration time
            tline = fgetl(fid);%#steps per measure
            tline = fgetl(fid);%#measures per dump
            tline = fgetl(fid);%#max dumps
            tline = fgetl(fid);%#blank
            tline = fgetl(fid);%#Spin Boson Parameters
            tline = fgetl(fid);%#g
            tline = fgetl(fid);%#alpha
            C = textscan(tline, '%f');
            obj.alpha = C{1}(1);
            
            tline = fgetl(fid);%#delta
            C = textscan(tline, '%f');
            obj.delta = C{1}(1);
            
            tline = fgetl(fid);%#v
            C = textscan(tline, '%f');
            obj.v = C{1}(1);
            
            tline = fgetl(fid);%#omega_c
            C = textscan(tline, '%f');
            obj.omega_c = C{1}(1);
            
            tline = fgetl(fid);%#blank
            tline = fgetl(fid);%#MPI/Parallel tempering params
            tline = fgetl(fid);%#parameter to vary
            C = textscan(tline, '%s');
            vary_param = C{1}{1};
            
            tline = fgetl(fid);%#variation in chosen param
            C = textscan(tline, '%f');
            param_delta = C{1}(1);
            
            %%Set q-c mapping and independent params
            obj.tc = obj.beta/obj.ly;
            obj.gamma = -0.5*log(obj.tc*obj.delta);
            obj.tcj = obj.tcj*obj.tc;
            obj.h = obj.h*obj.tc;
            
            if(strcmp('alpha', vary_param))
                obj.alpha = obj.alpha:param_delta:(obj.alpha + (obj.n_sims - 1)*param_delta);
            elseif(strcmp('r', vary_param))
                obj.a = obj.a:param_delta:(obj.a + (obj.n_sims - 1)*param_delta);
            end
        end
    end
    
end

