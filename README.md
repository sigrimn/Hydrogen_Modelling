# Flow Simulation and Hyperparameter Optimization for
Underground Hydrogen Storage

A 2D and a 3D model are implemented in the folder `workfiles`. All functions  
for optimization are either in the workfile scripts or the folder `objective
 functions`.  

Additional utility functions can be found in the `utils` folder.  
For a more detailed description of the implementations in this repository,
please refer to my specialization project.  
  

## How to Install and Run the Project

This project consists of MATLAB scripts. To set up the environment and run the code, follow these steps:

### 1️⃣ Clone this repository  
```sh
git clone https://github.com/sigrimn/Hydrogen_Modelling.git
cd Hydrogen_Modelling
```

### 2️⃣ Install and Clone Required Dependencies  

1. **Install MRST**  
   - Download and install MRST from [SINTEF's website](https://www.sintef.no/projectweb/mrst/download/).  

2. **Clone H2store Module**  
   - Clone the `H2store` module from Elyes Ahmed's repository:  
   ```sh
   git clone https://github.com/ElyesAhmed/MRST/tree/hydrogen/modules/H2store
   ```
   - Make sure this module is accessible within your MATLAB environment.

### 3️⃣ (Optional) Running on the MARKOV Server  

If you wish to run the project on the **MARKOV server** at the Department of Mathematical Sciences (NTNU), refer to the official guide:  
🔗 [MARKOV Server Wiki](https://wiki.math.ntnu.no/drift/stud/ommarkov/markov)  

#### **Setup on MARKOV**  
- Create a **work** location on the server where you store all **MRST and H2store functionalities**, as well as your output folder. Ask the techical group  
- Keep your own project files in a **home** location.  

#### **Advantages & Disadvantages**  
✅ **Advantage**: The `work` folder is directly on the MARKOV server, making execution faster.  
❌ **Disadvantage**: No backups are taken from the `work` location—only `home` is backed up.  


