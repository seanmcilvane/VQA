from qiskit import Aer, QuantumCircuit, ClassicalRegister, QuantumRegister, execute, assemble, transpile
import numpy as np
from qiskit.providers.aer import QasmSimulator
from qiskit.algorithms.optimizers import COBYLA
import math
from qiskit.visualization import plot_histogram
from utils import prob_dist


class VQA:

    def __init__(self, target_distribution, optimizer, backend, gate = "U3", layers = 4, entanglement = "Linear", shots = 1024):
        
        self.target_distr = target_distribution
        self.gate = gate
        self.layers = layers
        self.entanglement = entanglement
        self.optimizer = optimizer
        self.shots = shots
        self.backend = backend 
        self.n = int(math.log(len(target_distribution) ,2))
        
        
        if gate == "U3":
            m = np.random.rand(3*self.n*layers)
            self.params = m
            self.initial_params = m
        if gate == "RYRZ":
            m = np.random.rand(3*self.n*layers)
            self.params = m
            self.initial_params = m
        
        
    def Ansatz(self):
    
        qr = QuantumRegister(self.n)
        cr = ClassicalRegister(self.n)
        qc = QuantumCircuit(qr, cr)
        i = 0

 ########################################################################  
    
        if self.gate == "U3":
        
            for d in range(self.layers):
            
                for j in range(self.n):
                
                    qc.u(self.params[i+j], self.params[i+j+1], self.params[i+j+2], qr[j])
                i += self.n
                
                if d == self.layers - 1:
                    break
            
                if self.entanglement == "Linear":
            
                    for k in range(self.n):
                        if k == self.n-1:
                            break
                        qc.cx(qr[k], qr[k+1])
                
                elif self.entanglement == "Full":
            
                    for k in range(self.n):
                        if k == self.n-1:
                            break
                
                        for p in range(self.n):
                    
                        # CNOT gates only acting on qubits after the control
                            if p <= k:
                                continue
                        
                            qc.cx(qr[k], qr[p])    
                else:
                    print("Specify entanglement as Full or Linear")
            qc.measure(qr,cr)
            
 ########################################################################           
    
        if self.gate == "RYRZ":
    
            for d in range(self.layers):
        
                for j in range(self.n):
                    qc.ry(self.params[j+i],qr[j])
        
                i += self.n 
                for j in range(self.n):
                    qc.rz(self.params[j+i], qr[j])
        
                i += self.n
        
                if d == self.layers - 1:
                    break
                
 ########################################################################  
 ########################################################################  

                if self.entanglement == "Linear":
            
                    for k in range(self.n):
                        if k == self.n-1:
                            break
                        qc.cx(qr[k], qr[k+1])
                
 ########################################################################  
    
                elif self.entanglement == "Full":
            
                    for k in range(self.n):
                        if k == self.n-1:
                            break
                
                        for p in range(self.n):
                    
                    # CNOT gates only acting on qubits after the control
                            if p <= k:
                                continue
                        
                            qc.cx(qr[k], qr[p])    
                else:
                    print("Specify entanglement as Full or Linear")
            qc.measure(qr,cr)
    
        return qc
    
    
    def objective_function(self, params):
            
        qc = self.Ansatz()
        t_qc = transpile(qc, self.backend)
        qobj = assemble(t_qc, shots = self.shots)
        result = self.backend.run(qobj).result()
    
        # Sort counts by key for output distribution so the cost function is calculated correctly
        unsorted_output = result.get_counts(qc)
 
        # Adapt the output dictionary to be compatible
        if len(unsorted_output) != len(self.target_distr):
    
            for j in range(len(self.target_distr)):  # Run for all possible states
            
                binary = format(j, "b")
                check = str(binary)
                if len(check) < self.n: 
                    while len(check) < self.n:
                        check = "0" + check # 

                if check in unsorted_output:          # Checks to see if a state is observed 
                    continue                          # Continue if that state is observed
                else:
                    unsorted_output[check] = 0.0      # Create dictionary entry for unobserved states and assign a count value of 0
                            
    
    
        sorted_output = dict(sorted(unsorted_output.items()))
 
        #Convert obtained counts into probability vector
        output_distr = prob_dist(sorted_output, self.shots)
    
        # Calculate the cost the sum of the distances between the output distribution and the target distribution of each state 
    
        cost = sum([np.abs(output_distr[i] - self.target_distr[i]) for i in range(len(self.target_distr))])
    
        return cost
    
    def find_target_distribution(self):
        
        ret = self.optimizer.optimize(num_vars = len(self.params), objective_function = self.objective_function, initial_point = self.params)
        
        self.params = ret[0]
                                  

        qc = self.Ansatz()
        t_qc = transpile(qc, self.backend)
        qobj = assemble(t_qc, shots = self.shots)


        # get_counts returns an unsorted dictionary of the counts as values and the state as keys. This sorts
        # the dictionay by key so that the output distribution will be in the same order as the target_distribution.
        # The ordered dictionary is for example: 
        # {'0000', count of that state, '0001', count of that state, ... , '1111', count of that state}

        unsorted_output = self.backend.run(qobj).result().get_counts(qc)
        self.sorted_output = dict(sorted(unsorted_output.items()))
        output_distr = prob_dist(self.sorted_output, self.shots)


        print("Target Distribution:", self.target_distr)
        print("Obtained Distribution:", output_distr)
        print("Output Error (Manhattan Distance):", ret[1])
        print("Parameters Found:", ret[0])
        print("Intitial Parameters:", self.initial_params)
        print("Repetitions:", ret[2]) # Print how many repetitions the optimer needed to reach the tolerance
    
        self.output_distr = output_distr
        self.found_parameters = ret[0]
        
    
    
    
    # Create dictionary for target distribution with digital representation of states as keys and their respective expected counts
# as the corresponding values, since plot_histogram takes dictionary as input

# Creates keys for target dictionary that so that it matches with keys of output distribution (EX: 10 -> 010 for a 3 qubit system)

    def plot(self):
        
        number = []
        for i in range(len(self.target_distr)):
            binary = format(i, "b") # converts decimal to binary
            string = str(binary)
            if len(string) < self.n: 
                while len(string) < self.n:
                    string = "0" + string

            number.append(string)
    
        self.target_dictionary = dict(zip(number, self.shots*self.target_distr))

        # Plot histogram that compares target and output distributions
        figure = plot_histogram([self.target_dictionary, self.sorted_output], figsize = (20,8), legend = ['Target Dist', 'Output Dist'])


        return plot_histogram([self.target_dictionary, self.sorted_output], figsize = (20,8), legend = ['Target Dist', 'Output Dist'])
    
    
    