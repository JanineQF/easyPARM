#!/bin/bash
###################################################################################################################
#   Automated Force Fields for Metals     /$$$$$$$   /$$$$$$  /$$$$$$$  /$$      /$$                              # 
#                                        | $$__  $$ /$$__  $$| $$__  $$| $$$    /$$$                              #
#   /$$$$$$   /$$$$$$   /$$$$$$$ /$$   /$$| $$  \ $$| $$  \ $$| $$  \ $$| $$$$  /$$$$                             #
#  /$$__  $$ |____  $$ /$$_____/| $$  | $$| $$$$$$$/| $$$$$$$$| $$$$$$$/| $$ $$/$$ $$                             #
# | $$$$$$$$  /$$$$$$$|  $$$$$$ | $$  | $$| $$____/ | $$__  $$| $$__  $$| $$  $$$| $$                             #
# | $$_____/ /$$__  $$ \____  $$| $$  | $$| $$      | $$  | $$| $$  \ $$| $$\  $ | $$                             #
# |  $$$$$$$|  $$$$$$$ /$$$$$$$/|  $$$$$$$| $$      | $$  | $$| $$  | $$| $$ \/  | $$                             #
#  \_______/ \_______/|_______/  \____  $$|__/      |__/  |__/|__/  |__/|__/     |__/                             #
#                               /$$  | $$                                                                         #
#                              |  $$$$$$/              Ver. 4.25 - 19 May 2026                                    #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################


# Detect the directory where the current script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the directory from which easyPARM.sh was executed 
RUN_DIR="$1"

echo " "
echo "Running in directory: $RUN_DIR"
echo "Script directory: $SCRIPT_DIR"

# Ensure all outputs go to RUN_DIR
cd "$RUN_DIR"
# Remove files if it exists
files_to_remove=("dihedral.dat" "distance.dat" "esout" "atom_type.dat" "COMPLEX_modified.mol2" "COMPLEX_modified.frcmod" \
    "forcefield2.dat" "forcefield.dat" "metal_number.dat" "new_atomtype.dat" "temp_COMPLEX_modified.frcmod" \
    "temp.dat" "updated_COMPLEX_modified.frcmod" "updated_COMPLEX_modified2.frcmod" "angle.dat" "bond_angle_dihedral_data.dat" "new_atomtype1.dat"\
    "qout" "punch" "QOUT" "ATOMTYPE.INF" "leap.log" "updated_updated_COMPLEX_modified2.frcmod" "metals_complete.dat" "more_metal.dat" "new_atomtype2.dat" "REF_COMPLEX.mol2" "limited_data.dat" "line_number.dat" "ONE.mol2" "Reference_atom_type.dat" "COMPLEX.mol2" "ZEMA.mol2" "easyPARM_atomtype.dat" ) 

for file in "${files_to_remove[@]}"; do
    if [ -e "$file" ]; then
        rm -f "$file"
    fi
done

# Function to ask the user for input again if the information provided was incorrect
get_valid_input() {
    local prompt="$1"
    local valid_options="$2"
    local user_input
    while true; do
        read -p "$prompt" user_input
        # Trim any leading or trailing whitespace
        user_input=$(echo "$user_input" | xargs)
        
        # Check if user_input is within valid_options (by checking exact match)
        if [[ " $valid_options " =~ " $user_input " ]]; then
            echo "$user_input"
            return
        else
            echo "Invalid input. Please try again." >&2
        fi
    done
}

# Handle user's choice for Amber
echo " "
echo "================================="
echo "  Amber Configuration Menu "
echo "================================="
echo "Select your option:"
echo "1- Use currently loaded Amber"
echo "2- Specify Amber installation path"

while true; do
    choice=$(get_valid_input "Enter your choice: " "1 2")
    
    if [ "$choice" = "1" ]; then
        echo " "
        echo "Amber is assumed to be already loaded. Skipping sourcing."
        break
    elif [ "$choice" = "2" ]; then
        # Prompt the user for the path to Amber
        echo " "
        read -p "Please provide the path for Amber: " amber_path
        
        # Check if the provided path exists and the amber.sh file exists
        if [ -f "${amber_path}/amber.sh" ]; then
            # Source the Amber environment
            source "${amber_path}/amber.sh"
            
            # Check if sourcing was successful
            if [ $? -ne 0 ]; then
                echo "Failed to source Amber. Please check the path."
                continue
            fi
            echo "Amber environment sourced successfully."
            break
        else
            echo "The provided path is incorrect or amber.sh is missing. Please try again."
        fi
    fi
done

# Function to get user input
get_user_input() {
    # Get total charge
    while true; do
	echo " "
        read -p "Please provide the total charge: " charge_total
        if [[ "$charge_total" =~ ^-?[0-9]+$ ]]; then
            break
        else
	    echo " "
            echo "Invalid input. Please enter an integer for the total charge."
        fi
    done

    while true; do
	echo " "
        read -p "Please provide the total multiplicity: " multi_total
        if [[ "$multi_total" =~ ^-?[1-9]+$ ]]; then
            break
        else
	    echo " "
            echo "Invalid input. Please enter an integer for the total multiplicity."
        fi
    done
    # Ask the user to provide the optimized XYZ geometry file
    while true; do
	echo " "
        read -p "Please provide the optimized XYZ geometry file: " xyz_file
        if [[ ! "$xyz_file" =~ \.xyz$ ]]; then
            echo "Error: The file must have a .xyz extension. Please try again."
            continue
        fi
        if [ ! -f "$RUN_DIR/$xyz_file" ]; then
            echo "Error: Optimized structure is not found in $RUN_DIR. Please check the file name and try again."
            continue
        fi
        break
    done

    # Execute the Python script to calculate bond, angles and dihedrals
    if [ -f "$SCRIPT_DIR/02_get_bond_angle.py" ]; then
        python3 "$SCRIPT_DIR/02_get_bond_angle.py" "$RUN_DIR/$xyz_file"
    else
        echo "Script 02_get_bond_angle.py not found in $SCRIPT_DIR. Exiting."
        exit 1
    fi

    # Check if the script was successful
    if [ $? -ne 0 ]; then
        echo "Failed to execute 02_get_bond_angle.py. Exiting."
        exit 1
    fi

    echo " "
    echo "================================="
    echo "  Charge Method Selection Menu "
    echo "================================="
    echo "Select the charge calculation method:"
    echo "2- ORCA (ORCA Fit Charges)"
    echo "3- ORCA (RESP Charges)"
    charge_output=$(get_valid_input "Enter your choice: " "2 3")


    if [ "$charge_output" -eq 2 ]; then
        get_orca_input
    elif [ "$charge_output" -eq 3 ]; then
        get_orca_input
    fi
}


get_orca_input() {
    get_atom_type
    run_antechamber_orca
}


get_atom_type() {
    echo " "
    echo "================================="
    echo "  Atom Type Selection Menu "
    echo "================================="
    echo "Please select the atom type: "
    echo "1- Amber Force Field (AMBER)"
    echo "2- General Amber Force Field (GAFF)"
    echo "3- General Amber Force Field (GAFF2)"
    atom_type=$(get_valid_input "Enter your choice: " "1 2 3")

    # Map atom type to appropriate flag
    case $atom_type in
        1) at_type="amber";;
        2) at_type="gaff";;
        3) at_type="gaff2";;
    esac
}


run_antechamber_orca() {
    # Convert xyz to pdb
    if [ -f "$SCRIPT_DIR/xyz_to_pdb.py" ]; then
        python3 "$SCRIPT_DIR/xyz_to_pdb.py" "$RUN_DIR/$xyz_file" "$RUN_DIR/COMPLEX.pdb"
    else
        echo "Script xyz_to_pdb.py not found in $SCRIPT_DIR. Exiting."
        exit 1
    fi
    
    # Check if the script was successful
    if [ $? -ne 0 ]; then
        echo "Failed to execute xyz_to_pdb.py. Exiting."
        exit 1
    fi
    
    # Read charge output
    while true; do
        echo " "
        if [ "$charge_output" -eq 2 ]; then
		read -p "Please provide the charge output file (e.g., .log, .out): " charge_data
		if [[ ! "$charge_data" =~ \.(log|out)$ ]]; then
			echo "Error: The file must have a .log or .out extension. Please try again."
			continue
		fi
		if [ ! -f "$RUN_DIR/$charge_data" ]; then
			echo "Charge data file not found in $RUN_DIR. Please check the file name and try again."
			continue
		fi
	elif [ "$charge_output" -eq 3 ]; then
		read -p "Please provide the charge output file generated by weasel ( .txt): e.g. structure_charges.txt " charge_data
		if [[ ! "$charge_data" =~ \.txt$ ]]; then
			echo "Error: The file must have a .txt extension. Please try again."
			continue
		fi
		if [ ! -f "$RUN_DIR/$charge_data" ]; then
			echo "Charge data file not found in $RUN_DIR. Please check the file name and try again."
			continue
		fi

	fi
        break
    done
    
    # ATOMTYPE: AMBER
    # First antechamber attempt
    antechamber -i "$RUN_DIR/COMPLEX.pdb" -fi pdb -o "$RUN_DIR/COMPLEX.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -at "$at_type" -dr no -j 5 > "$RUN_DIR/temp.dat" 2>&1

    # ATOMTYPE GAFF
    # Apply our approach to detect the atom type.
    if [[ "$atom_type" -eq 2 || "$atom_type" -eq 3 ]]; then
	    
	    # Revise the mol2
	    python3 "$SCRIPT_DIR/03_correct_mol2.py" "$RUN_DIR"
	    # Update the bond section in mol2
	    python3 "$SCRIPT_DIR/atomtype_helper.py" "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/distance_type.dat" "$RUN_DIR/COMREF.mol2" > "$RUN_DIR/temp.dat" 2>&1
	    # Detect the atom type
	    python3 "$SCRIPT_DIR/atomtype_detector.py" "$RUN_DIR/COMREF.mol2" "$RUN_DIR/distance.dat" "$RUN_DIR/angle.dat" > "$RUN_DIR/temp.dat" 2>&1 
	    if [ -f "$RUN_DIR/easyPARM.mol2" ]; then
		    cp "$RUN_DIR/easyPARM.mol2" "$RUN_DIR/COMPLEX.mol2"
	    fi
    fi
    
    # Revise Atom Type
    python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
    
    # Additional processing steps: only necessary if number of metals larger than 5: warum auch immer ... 
    if [ ! -f "$RUN_DIR/limited_data.dat" ]; then
        if [ -f "$RUN_DIR/no_metal.dat" ]; then
            antechamber -i "$RUN_DIR/COMPLEX.pdb" -fi pdb -o "$RUN_DIR/COMPLEX.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -at "$at_type" > "$RUN_DIR/temp.dat" 2>&1
        else 
            python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
            antechamber -i "$RUN_DIR/mol.pdb" -fi pdb -o "$RUN_DIR/ONE.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -at "$at_type" -dr no > "$RUN_DIR/temp.dat" 2>&1
            python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
            if [ -f "$RUN_DIR/ONE.mol2" ]; then
                mv "$RUN_DIR/COMPLEX_modified.mol2" "$RUN_DIR/COMPLEX.mol2"
            fi
	    antechamber -i "$RUN_DIR/COMPLEX.pdb" -fi pdb -o "$RUN_DIR/ONE2.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -at "$at_type" -dr no > "$RUN_DIR/temp.dat" 2>&1
	    if [ -f "$RUN_DIR/ONE2.mol2" ]; then
                    cp "$RUN_DIR/ONE2.mol2" "$RUN_DIR/COMPLEX.mol2"
            fi
        fi
    else 
        break 
    fi
    
    # Check if COMPLEX.mol2 was generated after all attempts
    if [ ! -f "$RUN_DIR/COMPLEX.mol2" ]; then
        echo " "
	echo "Failed to generate COMPLEX.mol2 after multiple attempts."
	echo -e "\n\033[0;31mPlease verify the availability of Antechamber.\033[0m"
	echo " " 
        retry=$(get_valid_input "Would you like to provide a different charge output file? (y/n)" "y n")
        case "$retry" in
            [yY]) 
                echo "Retrying with a different charge output file..."
                run_antechamber_orca
                return
                ;;
            [nN]) 
                echo "Exiting due to failed antechamber command."
                exit 1
                ;;
        esac
    else
        echo "Antechamber command executed successfully. "
    fi

    # Update mol2 file with CHELPG charges (if selected)and Weasel RESP charges
    if [ "$charge_output" -eq 2 ] ; then
        awk '/ Charges/,/Total charge:/' "$RUN_DIR/$charge_data" | \
            grep -E '^\s*[0-9]+' | \
            awk '{print $NF}' > "$RUN_DIR/charges.dat"
        cp "$RUN_DIR/charges.dat" "$RUN_DIR/charges.chg"
       	python3 "$SCRIPT_DIR/Retrieve_RESP_Charges.py"
    elif [ "$charge_output" -eq 3 ]; then
        tail -n +2 "$RUN_DIR/$charge_data" | sed '/^[[:space:]]*$/d' > "charges.chg"
        python3 "$SCRIPT_DIR/Retrieve_RESP_Charges.py" 
    fi
}


# Main execution
get_user_input

# Ask the user for the input file
echo " "
echo "==========================================="
echo "  Frequencies Format Selection Menu "
echo "==========================================="
echo "Please select the format you will provide:"
echo "1- Orca Output"


qm_output=$(get_valid_input "Enter your choice: " "1")

case $qm_output in
    1)
        echo " "
	echo "You've selected Orca Output"
	if [ "$charge_output" -eq 1 -o "$charge_output" -eq 3 -o "$charge_output" -eq 4 -o "$charge_output" -eq 5 ]; then
            while true; do
                echo " "
                read -p "Please provide the Orca output file (.log or .out): " orca_output
                if [ ! -f "$RUN_DIR/$orca_output" ]; then
                    echo "Orca output file not found in $RUN_DIR. Please check the file name and try again."
                    continue
                fi
                break
            done
            while true; do
                echo " "
                read -p "Please provide the Orca hessian file (.hess): " orca_hessian
                if [ ! -f "$RUN_DIR/$orca_hessian" ]; then
                    echo "Orca hessian file not found in $RUN_DIR. Please check the file name and try again."
                    continue
                fi
                break
            done
            python3 "$SCRIPT_DIR/Seminario_method_ORCA.py" "$RUN_DIR/$orca_hessian" "$RUN_DIR/$orca_output" > "$RUN_DIR/temp.dat"
	    
	    if [ $? -ne 0 ]; then
		    echo "Failed to execute Seminario_method_ORCA.py. Exiting."
		    echo -e "\n\033[0;31mPlease check the output.\033[0m"

		    exit 1
	    fi

        elif [ "$charge_output" -eq 2 ]; then
            while true; do
                echo " "
                read -p "Please provide the Orca hessian file (.hess): " orca_hessian
                if [ ! -f "$RUN_DIR/$orca_hessian" ]; then
                    echo "Orca hessian file not found in $RUN_DIR. Please check the file name and try again."
                    continue
                fi
                break
            done
            python3 "$SCRIPT_DIR/Seminario_method_ORCA.py" "$RUN_DIR/$orca_hessian" "$RUN_DIR/$charge_data" > "$RUN_DIR/temp.dat"
	    if [ $? -ne 0 ]; then
		    echo "Failed to execute Seminario_method_ORCA.py. Exiting."
		    echo -e "\n\033Please verify the following:\033[0m"
		    echo -e "  1. Check if the calculation terminated correctly"
		    echo -e "  2. Verify frequency calculation completed separately"

		    exit 1
	    fi
        fi
        ;;
    
esac
    

#Remove unnecessary file
if [  -f "$RUN_DIR/metals_complete.dat" ]; then
    rm "$RUN_DIR/metals_complete.dat" 
fi

# Run the Python and shell scripts in sequence with checks
# Generate a .mol2 file
# Correct the .mol2 file
# Generate a .frcmod file
# Assign new atom types in the .mol2 and .frcmod file
#for script in 03_correct_mol2.py correct_atom_type.py 04_parmch2_frcmod.sh; do
for script in 03_correct_mol2.py correct_atom_type.py 04_parmch2_frcmod.sh; do

    if [[ $script == *.py ]]; then
        if [ -f "$SCRIPT_DIR/$script" ]; then
            python3 "$SCRIPT_DIR/$script" "$RUN_DIR"
	    if [  -f "$RUN_DIR/COMPLEX_modified_atom_type.mol2" ]; then
           # mv "$RUN_DIR/COMPLEX_modified_atom_type.mol2" "$RUN_DIR/COMPLEX.mol2" 
		    cp "$RUN_DIR/COMPLEX_modified_atom_type.mol2" "$RUN_DIR/COMPLEX.mol2"
	    fi

            if [ $? -ne 0 ]; then
                echo "Failed to execute $script. Exiting."
                exit 1
            fi
        else
            echo "Script $script not found in $SCRIPT_DIR. Exiting."
            exit 1
        fi
    elif [[ $script == *.sh ]]; then
        if [ -f "$SCRIPT_DIR/$script" ]; then
            bash "$SCRIPT_DIR/$script" "$RUN_DIR" "$atom_type"
            if [ $? -ne 0 ]; then
                echo "Failed to execute $script. Exiting."
                exit 1
            fi
        else
            echo "Script $script not found in $SCRIPT_DIR. Exiting."
            exit 1
        fi
    fi
done

echo "Done creating COMPLEX.mol2 file with assigned atom types and charges."


# Function to process metalloprotein with selected force field
process_metalloprotein() {
    local protein_ff=$1

    if [ ! -f "$SCRIPT_DIR/metalloprotein.py" ]; then
        echo "Error: metalloprotein.py not found in $SCRIPT_DIR" >&2
        return 1
    fi

    case "$protein_ff" in
        1) lib_file="fb15.lib" ;;
        2) lib_file="ff12SB.lib" ;;
        3) lib_file="ff14SB.lib" ;;
        4) lib_file="ff19SB.lib" ;;
        *)
            echo "Error: Invalid protein force field selection: $protein_ff" >&2
            echo "Valid options are 1-4" >&2
            return 1
            ;;
    esac

    python3 "$SCRIPT_DIR/metalloprotein.py" "$SCRIPT_DIR/libraries/$lib_file"
    mv "$RUN_DIR/COMPLEX_updated.mol2" "$RUN_DIR/easyCOMPLEX.mol2"
}

metalloprotein_choice=$(get_valid_input "Does your structure belong to MetalloProtein ? (y/n): " "y n yes no Y N YES NO Yes No")
echo " "

if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then

    echo " "
    while true; do
    
    	cp "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/NEW_COMPLEX.mol2" 
    	if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then
		
		read -p "Please provide the metalloprotein pdb file: " protein_pdb
	        
		echo " "
		echo "======================================================"
		echo "        Protein Force Field Integration Menu"
		echo "======================================================"
		echo "Select the protein FF to be used for the metal-coordinated standard amino-acid residues"
		echo " "
		echo "1- fb15"
		echo "2- ff12SB "
		echo "3- ff14SB"
		echo "4- ff19SB"
		echo " "
		protein_ff=$(get_valid_input "Enter your choice: " "1 2 3 4")

		if [ ! -f "$RUN_DIR/$protein_pdb" ]; then
		    echo "Metalloprotein file not found in $RUN_DIR. Please check the file name and try again."
		    continue
		fi

		pdb4amber -i "$RUN_DIR/$protein_pdb" -o "$RUN_DIR/metalloprotein_easyPARM.pdb"  > "$RUN_DIR/temp.dat" 2>&1 
		if [ -f "$SCRIPT_DIR/metalloprotein.py" ]; then
			process_metalloprotein "$protein_ff"
		else
			echo "Script metalloprotein.py not found in $SCRIPT_DIR. Exiting."
			exit 1
		fi

		# Check if the script was successful
		if [ $? -ne 0 ]; then
			echo "Failed to execute metalloprotein.py. Exiting."
			exit 1
		fi
    	fi
    	
	
	mv "$RUN_DIR/easyCOMPLEX.mol2" "$RUN_DIR/COMPLEX_modified.mol2" 
    	
    	if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then
	python3 "$SCRIPT_DIR/distribute_metalloprotein_charge.py" "$RUN_DIR/$protein_pdb" "$RUN_DIR/COMPLEX_modified.mol2" $charge_total 
	fi 
    	

    	
	mv "$RUN_DIR/updated_easy_COMPLEX.mol2" "$RUN_DIR/COMPLEX_modified.mol2" 
    	break
    done
else
	if [ -f "$RUN_DIR/COMPLEX.lib" ]; then
		rm "$RUN_DIR/COMPLEX.lib"
    	fi
    
fi



# Select appropriate preparation script in order to generate the frcmod depend on the complexation of the system
# more_metal.dat refer that there are more than metal or non standard atom then the code will select the code depend on that
if [ -f "$RUN_DIR/more_metal.dat" ]; then

    python3 "$SCRIPT_DIR/05_prepare_mol2_frcmod_more_atom.py" "$RUN_DIR"
    cp "$RUN_DIR/new_atomtype.dat" "$RUN_DIR/new_atomtype1.dat"    
    cp "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/REF_COMPLEX.mol2"    
    cp "$RUN_DIR/NEW_COMPLEX.mol2" "$RUN_DIR/COMPLEX.mol2"    
    mv "$RUN_DIR/metal_number.dat" "$RUN_DIR/metals_complete.dat"
    rm "$RUN_DIR/more_metal.dat"
    python3 "$SCRIPT_DIR/generate_preforcefield.py" "$RUN_DIR"
    if [ $? -ne 0 ]; then
        echo "Failed to execute 05_prepare_mol2_frcmod_more_atom.py. Exiting."
        exit 1
    fi
fi

if [ ! -f "$RUN_DIR/more_metal.dat" ]; then
    python3 "$SCRIPT_DIR/05_prepare_mol2_frcmod.py" "$RUN_DIR"
    if [ $? -ne 0 ]; then
        echo "Failed to execute 05_prepare_mol2_frcmod.py. Exiting."
        exit 1
    fi
fi

if [  -f "$RUN_DIR/metals_complete.dat" ]; then
    cp "$RUN_DIR/metals_complete.dat" "$RUN_DIR/metal_number.dat"
fi

if [  -f "$RUN_DIR/new_atomtype1.dat" ]; then
    cp "$RUN_DIR/new_atomtype.dat" "$RUN_DIR/new_atomtype2.dat"    
    cp "$RUN_DIR/REF_COMPLEX.mol2" "$RUN_DIR/COMPLEX.mol2"
    cat "$RUN_DIR/new_atomtype1.dat" >> "$RUN_DIR/new_atomtype2.dat"
fi

echo " "
if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then

	cp "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/1COMPLEX.mol2"
	if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then
		if [ -f "$SCRIPT_DIR/metalloprotein.py" ]; then
			process_metalloprotein "$protein_ff"
		else
			echo "Script metalloprotein.py not found in $SCRIPT_DIR. Exiting."
			exit 1
		fi

		# Check if the script was successful
		if [ $? -ne 0 ]; then
			echo "Failed to execute metalloprotein.py. Exiting."
			exit 1
		fi
	fi


	python3 "$SCRIPT_DIR/xyz_to_pdb.py" "$RUN_DIR/part_QM.xyz" "$RUN_DIR/part_QM.pdb"
	python3 "$SCRIPT_DIR/xyz_to_pdb.py" "$RUN_DIR/qm.xyz" "$RUN_DIR/qm.pdb"

	antechamber -i "$RUN_DIR/qm.pdb" -fi pdb -o "$RUN_DIR/QM.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -at "$at_type" -dr no -j 5 > "$RUN_DIR/temp.dat" 2>&1
	cp "$RUN_DIR/COMPLEX.pdb" "$RUN_DIR/backupCOMPLEX.pdb"
	cp "$RUN_DIR/qm.pdb" "$RUN_DIR/COMPLEX.pdb"
	cp "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/backupCOMPLEX.mol2"
	cp "$RUN_DIR/QM.mol2" "$RUN_DIR/COMPLEX.mol2"
	
	if [ -f "$RUN_DIR/limited_data.dat" ]; then
		rm -f "$RUN_DIR/limited_data.dat" 
	fi
    
	if [ -f "$RUN_DIR/line_number.dat" ]; then
		rm -f "$RUN_DIR/line_number.dat" 
	fi
	
	python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
	if [ ! -f "$RUN_DIR/limited_data.dat" ]; then
            python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
            antechamber -i "$RUN_DIR/mol.pdb" -fi pdb -o "$RUN_DIR/ONE.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -at "$at_type" -dr no > "$RUN_DIR/temp.dat" 2>&1
            python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
            if [ -f "$RUN_DIR/ONE.mol2" ]; then
                mv "$RUN_DIR/COMPLEX_modified.mol2" "$RUN_DIR/QM.mol2"
            fi
	    antechamber -i "$RUN_DIR/COMPLEX.pdb" -fi pdb -o "$RUN_DIR/ONE2.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -at "$at_type" -dr no > "$RUN_DIR/temp.dat" 2>&1
	    if [ -f "$RUN_DIR/ONE2.mol2" ]; then
		    cp "$RUN_DIR/ONE2.mol2" "$RUN_DIR/QM.mol2"
            fi
    	else 
        	break 
    	fi
    
       	mv "$RUN_DIR/backupCOMPLEX.pdb" "$RUN_DIR/COMPLEX.pdb"
	mv "$RUN_DIR/backupCOMPLEX.mol2" "$RUN_DIR/COMPLEX.mol2"
    
    	sed -i'' '$d' "$RUN_DIR/nonstand.pdb"
    	cat "$RUN_DIR/part_QM.pdb" >> "$RUN_DIR/nonstand.pdb"
    	cat "$RUN_DIR/part_QM.pdb" >> "$RUN_DIR/easynonstands.pdb"
    	pdb4amber -i "$RUN_DIR/easynonstands.pdb" -o "$RUN_DIR/easyPARM.pdb"  > "$RUN_DIR/temp.dat" 2>&1 
     
    	input_file="$RUN_DIR/easyPARM_residues.dat"

    	# Read each line from the input file
    	while IFS=' ' read -r pdbout mol2out residue_name; do
    	# Skip empty lines or lines starting with #
    		if [[ -z "$pdbout" || "$pdbout" == \#* ]]; then
        		continue
    		fi
    
    	# Run the antechamber command
    		antechamber -i "$pdbout" -fi pdb -o "$mol2out" -fo mol2 -s 2 -rn "$residue_name"  
   	done < "$input_file" > "$RUN_DIR/temp.dat" 2>&1

     
    	if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then
		python3 "$SCRIPT_DIR/update_metalloprotein_charge.py" > "$RUN_DIR/temp.dat" 2>&1
	fi 
    	

    	cp "$RUN_DIR/QM.mol2" "$RUN_DIR/COMPLEX.mol2"
    	python3 "$SCRIPT_DIR/02_get_bond_angle.py" "$RUN_DIR/qm.xyz"
    	python3 "$SCRIPT_DIR/03_correct_mol2.py" 
    	mv "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/QM.mol2"
    	mv "$RUN_DIR/1COMPLEX.mol2" "$RUN_DIR/COMPLEX.mol2"
    	mv "$RUN_DIR/easyCOMPLEX.mol2" "$RUN_DIR/COMPLEX_modified.mol2" 
    	python3 "$SCRIPT_DIR/02_get_bond_angle.py" "$RUN_DIR/$xyz_file"
else
    :
fi


# Generate the frcmod file with all the correct info and with a clean version
for script in 06_get_atom_type.py 07_Seminario_forcefield.py 08_update_forcefield.py 09_clean_updatedforcefield.py 10_postclean_updatedforcefield.py 11_retrieve_uffdata.py 13_final_clean.py; do
    if [ -f "$SCRIPT_DIR/$script" ]; then
        python3 "$SCRIPT_DIR/$script" "$RUN_DIR"

        if [ $? -ne 0 ]; then
            echo "Failed to execute $script. Exiting."
            exit 1
        fi
    else
        echo "Script $script not found in $SCRIPT_DIR. Exiting."
        exit 1
    fi
done


# Change the name of output
mv spaced_filtered_COMPLEX_modified2.frcmod COMPLEX.frcmod
cp NEW_COMPLEX.mol2 COMPLEX.mol2 

python3 "$SCRIPT_DIR/update_if_needed.py" "$RUN_DIR/COMPLEX.frcmod"

# Function to process metalloprotein parameters (frcmod)
process_metalloprotein_parm() {
    local protein_ff=$1
    
    if [ ! -f "$SCRIPT_DIR/metalloprotein_parm.py" ]; then
        echo "Error: metalloprotein_parm.py not found in $SCRIPT_DIR" >&2
        return 1
    fi
    
    case "$protein_ff" in
        1) frcmod_file="fb15.frcmod" ;;
        2) frcmod_file="ff12SB.frcmod" ;;
        3) frcmod_file="ff14SB.frcmod" ;;
        4) frcmod_file="ff19SB.frcmod" ;;
        *) 
            echo "Error: Invalid protein force field selection: $protein_ff" >&2
            echo "Valid options are 1-4" >&2
            return 1
            ;;
    esac
    
    python3 "$SCRIPT_DIR/metalloprotein_parm.py" "$SCRIPT_DIR/libraries/$frcmod_file" > "$RUN_DIR/temp.dat" 2>&1
}



resid_ID=$(get_valid_input "Would you like to change the residue ID (Default= mol)? (y/n): " "y n yes no Y N YES NO Yes No")
if [[ "${resid_ID,,}" =~ ^(y|yes)$ ]]; then
    read -r -p "Please provide the residue name: " resid_name
fi

if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then

    
    if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then
	    process_metalloprotein_parm "$protein_ff" 
	    mv "$RUN_DIR/nonstand.pdb" "$RUN_DIR/easyPARM_MetalloProtein.pdb"
	    pdb4amber -i "$RUN_DIR/easyPARM_MetalloProtein.pdb" -o "$RUN_DIR/easyPARM_MetalloProtein2.pdb"  > "$RUN_DIR/temp.dat" 2>&1 
	    mv "$RUN_DIR/easyPARM_MetalloProtein2.pdb" "$RUN_DIR/easyPARM_MetalloProtein.pdb" 
    fi

    
    mv "$RUN_DIR/QM.mol2" "$RUN_DIR/METAL.mol2"
    mv "$RUN_DIR/coordinated_residues.txt" "$RUN_DIR/Bond_Info.dat"
    echo "  "
    

    echo "=========================================="
    echo "         	Output Files "
    echo "=========================================="

    input_file="$RUN_DIR/easyPARM_residues.dat"
    line_number=1

# Read each line from the input file
    while IFS=' ' read -r pdbout mol2out residue_name; do
    	# Skip empty lines or lines starting with #
    	if [[ -z "$mol2out" || "$pdbout" == \#* ]]; then
        	continue
    	fi	
    
    # Print the mol2out with the line number
    	echo "Mol2 $line_number			  : $mol2out"
    	((line_number++))
    done < "$input_file" 
    if [[ "${resid_ID,,}" =~ ^(y|yes)$ ]]; then
	    cp "$RUN_DIR/COMPLEX.frcmod" "$RUN_DIR/COMPLEX_${resid_name}.frcmod"
	    #mv "$RUN_DIR/COMPLEX.frcmod" "$RUN_DIR/easyPARM.frcmod"
	    mv "$RUN_DIR/Bond_Info.dat" "$RUN_DIR/Bond_Info_${resid_name}.dat" 
	    mv "$RUN_DIR/Hybridization_Info.dat" "$RUN_DIR/Hybridization_Info_${resid_name}.dat"
	    if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then
		    mv "$RUN_DIR/easyPARM_MetalloProtein.pdb" "$RUN_DIR/easyPARM_MetalloProtein_${resid_name}.pdb"
		    sed -i'' "s/\<mol\>/${resid_name}/g" "$RUN_DIR/easyPARM_MetalloProtein_${resid_name}.pdb"
	    fi

	    cp "$RUN_DIR/METAL.mol2" "$RUN_DIR/${resid_name}.mol2"
	    sed -i'' "s/\<mol\>/${resid_name}/g" "$RUN_DIR/${resid_name}.mol2"

	    echo "Mol2  		    	  : ${resid_name}.mol2"
	    echo "Frcmod                    : COMPLEX_${resid_name}.frcmod"
	    echo "Bond Information          : Bond_Info_${resid_name}.dat"
	    echo "Lib    		    	  : COMPLEX.lib"
	    echo "New Atom Type             : Hybridization_Info_${resid_name}.dat"
	    echo "Tleap Input   		  : Tleap.input"
	    if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then
		    echo "MetalloProtein pdb        : easyPARM_MetalloProtein_${resid_name}.pdb"
	    fi
    else
	    echo "Mol2                      : METAL.mol2"
	    echo "Bond Information          : Bond_Info.dat"
	    echo "Frcmod                    : COMPLEX.frcmod"
	    echo "New Atom Type             : Hybridization_Info.dat"
	    echo "Lib    		    	  : COMPLEX.lib"
	    echo "MetalloProtein pdb        : easyPARM_MetalloProtein.pdb"
	    echo "Tleap Input               : Tleap.input"
    fi
else 
    	if [[ "${resid_ID,,}" =~ ^(y|yes)$ ]]; then
	    #sed -i  -E "s/ mol /${resid_name}/g" "$RUN_DIR/COMPLEX.mol2"
	    	sed -i'' "s/\<mol\>/${resid_name}/g" "$RUN_DIR/COMPLEX.mol2"
	    	sed -i'' "s/\<mol\>/${resid_name}/g" "$RUN_DIR/COMPLEX.pdb"
    	fi

# Print all the output name
	echo "================================="
	echo "         Output Files "
	echo "================================="
	echo "Mol2   		: COMPLEX.mol2"
	echo "Frcmod 		: COMPLEX.frcmod"
	echo "PDB    		: COMPLEX.pdb"
	echo "Lib    		: COMPLEX.lib"
    	echo "New Atom Type	: Hybridization_Info.dat"
	echo "Tleap Input   	: tleap.input"

fi	

if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then

        # run for only the small residue 
	python3 "$SCRIPT_DIR/02_get_bond_angle.py" "$RUN_DIR/part_QM.xyz"

	# Prepare library input for tleap
	echo "source leaprc.gaff" > input_library.tleap
	echo "loadamberparams COMPLEX.frcmod" >> input_library.tleap
	echo "mol = loadmol2 METAL.mol2" >> input_library.tleap
	echo "check mol" >> input_library.tleap
	echo "charge mol" >> input_library.tleap
	echo "savepdb mol COMPLEX.pdb" >> input_library.tleap
	echo "saveoff mol COMPLEX.lib" >> input_library.tleap
	echo "quit" >> input_library.tleap
	tleap -f input_library.tleap > leap.log
	tleap -f ALL_RESIDUE_tleap.input > leap.log

	#Function to correct lib file
	python3 "$SCRIPT_DIR/12_generate_lib.py" > "$RUN_DIR/temp.dat" 2>&1 
	
	if [[ "${resid_ID,,}" =~ ^(y|yes)$ ]]; then 
		sed -i'' "s/\<mol\>/${resid_name}/g" "$RUN_DIR/COMPLEX.lib"
		rm "$RUN_DIR/COMPLEX.frcmod"
		rm "$RUN_DIR/METAL.mol2"
	fi
	# Generation of tleap input
	echo "source leaprc.gaff2" > tleap.input
	echo "source leaprc.water.tip3p" >> tleap.input
        echo " " >> tleap.input
	echo "loadamberparams frcmod.ionsjc_tip3p" >> tleap.input
        echo " " >> tleap.input
	echo "# Load force field for residues and metal" >> tleap.input
	if [[ "${resid_ID,,}" =~ ^(y|yes)$ ]]; then 
		echo "loadamberparams COMPLEX_${resid_name}.frcmod" >> tleap.input
	else
		echo "loadamberparams COMPLEX.frcmod" >> tleap.input
	fi
        echo " " >> tleap.input
	echo "# Modify the system force field to your preferred version; for example, use ff14SB instead of ff19SB." >> tleap.input
	if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]] ; then
		echo "source leaprc.protein.ff19SB " >> tleap.input
	fi

        echo " " >> tleap.input
	echo "# Load library for residues and metal" >> tleap.input
	echo "loadoff COMPLEX.lib" >> tleap.input 
	if [[ "${resid_ID,,}" =~ ^(y|yes)$ ]]; then 
        	echo " " >> tleap.input
		echo "PRO = loadpdb "easyPARM_MetalloProtein_${resid_name}.pdb"" >> tleap.input
        	echo " " >> tleap.input
		cat "$RUN_DIR/Bond_Info_${resid_name}.dat" >> tleap.input  
        	echo -e "\n" >> tleap.input
	else
        	echo " " >> tleap.input
		echo "PRO = loadpdb "easyPARM_MetalloProtein.pdb"" >> tleap.input
        	echo " " >> tleap.input
		cat "$RUN_DIR/Bond_Info.dat" >> tleap.input
        	echo -e "\n" >> tleap.input
	fi
        echo "# Save vacuum files" >> tleap.input
        echo "savepdb PRO System_vacuum.pdb" >> tleap.input
	echo "saveamberparm PRO System_vacuum.prmtop System_vacuum.inpcrd" >> tleap.input	
        echo " " >> tleap.input
	echo "# Solvate the system" >> tleap.input
	echo "solvateoct PRO TIP3PBOX 10.0" >> tleap.input
	echo "addions PRO Cl- 0. " >> tleap.input
	echo "addions PRO Na+ 0. " >> tleap.input
        echo " " >> tleap.input
	echo "# Save solvated files" >> tleap.input
	echo "saveamberparm PRO System_solvated.prmtop System_solvated.inpcrd" >> tleap.input
	echo "savepdb PRO System_solvated.pdb" >> tleap.input
	echo "quit" >> tleap.input
else 
	echo "source leaprc.gaff" > input_library.tleap
	echo "loadamberparams COMPLEX.frcmod" >> input_library.tleap
	echo "mol = loadmol2 "COMPLEX.mol2"" >> input_library.tleap
	echo "check mol" >> input_library.tleap
	echo "charge mol" >> input_library.tleap
	echo "savepdb mol COMPLEX.pdb" >> input_library.tleap
	echo "saveoff mol COMPLEX.lib" >> input_library.tleap
	echo "quit" >> input_library.tleap

	#Function to correct lib file
	tleap -f input_library.tleap > leap.log

	python3 "$SCRIPT_DIR/12_generate_lib.py" > "$RUN_DIR/temp.dat" 2>&1
	if [[ "${resid_ID,,}" =~ ^(y|yes)$ ]]; then 
		sed -i'' "s/\<mol\>/${resid_name}/g" "$RUN_DIR/COMPLEX.lib"
	fi
	# Generation of tleap input
	echo "source leaprc.gaff2" > tleap.input
	echo "source leaprc.water.tip3p" >> tleap.input
        echo " " >> tleap.input
	echo "loadamberparams frcmod.ionsjc_tip3p" >> tleap.input
        echo " " >> tleap.input
	echo "loadamberparams COMPLEX.frcmod" >> tleap.input
        echo " " >> tleap.input
	echo "# Uncomment the correct force field with respect to your system" >> tleap.input
	echo "# Modify the system force field to your preferred version; for example, use ff14SB instead of ff19SB." >> tleap.input
	echo "#source leaprc.protein.ff19SB " >> tleap.input
	echo "#source leaprc.DNA.OL24 " >> tleap.input
        echo " " >> tleap.input
	echo "# Load library for metal complex" >> tleap.input
	echo "loadoff COMPLEX.lib" >> tleap.input 
        echo " " >> tleap.input
	echo "PRO = loadpdb "Your_Whole_System.pdb"" >> tleap.input
        echo " " >> tleap.input
        echo "# Save vacuum files" >> tleap.input
        echo "savepdb PRO System_vacuum.pdb" >> tleap.input
	echo "saveamberparm PRO System_vacuum.prmtop System_vacuum.inpcrd" >> tleap.input	
        echo " " >> tleap.input
	echo "# Solvate the system" >> tleap.input
	echo "solvateoct PRO TIP3PBOX 10.0" >> tleap.input
	echo "addions PRO Cl- 0. " >> tleap.input
	echo "addions PRO Na+ 0. " >> tleap.input
        echo " " >> tleap.input
	echo "# Save solvated files" >> tleap.input
	echo "saveamberparm PRO System_solvated.prmtop System_solvated.inpcrd" >> tleap.input
	echo "savepdb PRO System_solvated.pdb" >> tleap.input
	echo "quit" >> tleap.input

fi



if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then
	rm -f COMPLEX.mol2
	rm -f COMPLEX.pdb	
	rm -f metalloprotein_easyPARM_*
	rm -f *_renum.txt
	rm -f *_sslink
	rm -f *_nonprot.pdb
	rm -f charge_*.dat
	rm -f recalculated_charges.dat
	rm -f processed_charges.dat
	rm -f original_charges.dat
	rm -f coordination_analysis.txt
	rm -f charge_statistics.txt	
	rm -f fixed_charges.dat
		
fi
# Remove the unnecessary file
rm ANTECHAMBER*
files_to_remove=("dihedral.dat" "distance.dat" "esout" "atom_type.dat" "COMPLEX_modified.mol2" "COMPLEX_modified.frcmod" \
    "complex.fchk" "forcefield2.dat" "metal_number.dat" "temp_COMPLEX_modified.frcmod" "new_atomtype.dat" \
    "temp.dat" "updated_COMPLEX_modified.frcmod" "updated_COMPLEX_modified2.frcmod" "angle.dat" "new_atomtype1.dat"\
    "qout" "punch" "QOUT" "ATOMTYPE.INF" "leap.log" "updated_updated_COMPLEX_modified2.frcmod" "metals_complete.dat" "more_metal.dat" "new_atomtype2.dat" "REF_COMPLEX.mol2" "limited_data.dat" "mol.pdb" "line_number.dat" "ONE.mol2" "Reference_atom_type.dat" "REFQM.pdb" "NEW_COMPLEX.mol2"\
    "QM.pdb" "nonstand.pdb" "easyPARM.pdb" "easynonstands.pdb" "part_QM.xyz" "part_QM.pdb" "charge_qm.dat" "metalloprotein.pdb" "metalloprotein_easyPARM.pdb" "charges_all.dat"  "easyPARM_residues.dat" "reference_structure.xyz" "ALL_RESIDUE_tleap.input" "ZEMA.mol2" "easyPARM_atomtype.dat" "easyPARM.mol2" "COMREF.mol2" "COMPLEX.frcmod.bak" "metalloprotein_atomtype.dat" "filtered_COMPLEX_modified2.frcmod" "capping_link_atoms.dat") 

for file in "${files_to_remove[@]}"; do
    if [ -e "$file" ]; then
        rm -f "$file"
    fi
done

# Function to check and report low force constants
check_force_constants() {
    local frcmod_file="$1"
    local low_constants=""
    local in_bond=0

    while IFS= read -r line; do
        if [[ "$line" =~ ^BOND ]]; then
            in_bond=1
            continue
        elif [[ "$line" =~ ^ANGLE ]]; then
            break  # stop at ANGLE section
        fi

        if (( in_bond )) && [[ "$line" =~ ^[A-Za-z] ]]; then
            # Capture full atom-pair and force constant from the whole line
            if [[ "$line" =~ ^([A-Za-z0-9*\'-]+([[:space:]]*-[[:space:]]*[A-Za-z0-9*\'-]+)+)[[:space:]]+([-0-9.]+)[[:space:]]+([-0-9.]+) ]]; then
                atom_pair="${BASH_REMATCH[1]}"
                force="${BASH_REMATCH[3]}"

                # Use bc for float comparison
                if (( $(echo "$force < 20" | bc -l) )); then
                    low_constants+="  ${atom_pair}  has force constant  ${force}"$'\n'
                fi
            fi
        fi
    done < "$frcmod_file"

    if [ -n "$low_constants" ]; then
        echo "========================================================================="
        echo -e "\n\033[0;31mAbnormally low bond force constants detected!\033[0m"
        echo "========================================================================="
        echo "$low_constants"
        echo -e "\n\033[0;31mThis indicates severely stretched bonds in your structure.\033[0m"
        echo -e "\n\033[0;31mCheck your geometry optimization and level of theory before proceeding.\033[0m"
        echo "========================================================================="
        echo ""
    fi
}

# Main logic
if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then
    if [[ "${resid_ID,,}" =~ ^(y|yes)$ ]]; then
        frcmod_file="$RUN_DIR/COMPLEX_${resid_name}.frcmod"
    else
        frcmod_file="$RUN_DIR/COMPLEX.frcmod"
    fi
else
    frcmod_file="$RUN_DIR/COMPLEX.frcmod"
fi

check_force_constants "$frcmod_file"


echo " "
echo -n "🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉"
echo -e "\n\033[1;36m🏁 Crossed the finish line successfully! 🏁\033[0m"
echo "🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉🎉"
