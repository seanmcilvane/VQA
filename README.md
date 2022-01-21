Sean Mcilvane
seanmcilvane@gmail.com 
sm892@duke.edu

# VQA_heuristic_approach
A Variational Quantum Algorithm that uses a heuristic (or hardware-efficient ansatz) approach to create a desired probability distribution from a random probability or pre-set probability distribution.

The program creates a hardware-efficient variational form based on inputs from the user. The user can select how many qubits to put into the variational form. The user is able to customize the variational form with different types of entanglement between qubits and the types of gates used. Currently, the user can choose between Linear and Full entanglement and the user can also choose between using RY RZ gates or U3 gates on each qubit. A layer is one set of gates with one set entangling CNOT gates. The user can also choose how many times a layer is repeated for their circuit. It should be noted, however, that increasing the layers of the circuit will increase the number of parameters, runtime, and error if ran on a real quantum computer. If too many parameters are used, the user might overfit their training and thus get worse results than if they were to use less paramters.
