from django.http import HttpResponse
from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from django.core.files.storage import FileSystemStorage
from pymongo import MongoClient
from django.conf import settings
import datetime
from django.contrib.auth import authenticate, login, logout
from django.shortcuts import render, redirect
from django.contrib import messages
from django.contrib.auth.hashers import make_password, check_password
from django.conf import settings
from django.http import HttpResponse
import pymongo
from django.contrib.auth.decorators import login_required



def aboutUS(request):
    return HttpResponse("<b>Welcome to AffPRo</b>")

def contactDetails(request, id):
    return HttpResponse(id)

def homePage(request):
    return render(request,"index.html")



#from rdkit import Chem
#from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import io
import base64
from rdkit.Chem import Descriptors

def generate_ligand_descriptors(smiles):
    # Convert the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        raise ValueError("Invalid SMILES string")
    
    # Calculate the required descriptors
    descriptors = {
        "MolecularWeight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
        "NumRings": Descriptors.RingCount(mol),
        #"VanDerWaalsVolume": Descriptors.VanDerWaalsVolume(mol),  # Approximation for Van der Waals Volume
        #"BalabanIndex": Descriptors.BalabanIndex(mol),
    }

    return descriptors

from rdkit.Chem import Draw
def visualize_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol)
            buffer = io.BytesIO()
            img.save(buffer, format="PNG")
            img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
            buffer.close()
            return f'<img src="data:image/png;base64,{img_str}" alt="SMILES Visualization"/>'
        else:
            return "<p>Invalid SMILES string</p>"
    except Exception as e:
        return f"<p>Error generating visualization: {e}</p>"




def visualize_3D_smiles(smiles):
    try:
        # Convert SMILES to molecule
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)  # Add hydrogens
        AllChem.EmbedMolecule(mol)  # Generate 3D coordinates
        AllChem.UFFOptimizeMolecule(mol)  # Optimize structure
        pdb_block = Chem.MolToPDBBlock(mol)  # Convert to PDB format

        # Generate 3D visualization using py3Dmol
        view = py3Dmol.view(width=700, height=300)
        view.addModel(pdb_block, "pdb")  # Add PDB model
        view.setStyle({'stick': {}})  # Use stick model
        view.zoomTo()  # Zoom to fit the molecule

        # Generate HTML for embedding in template
        return view._make_html()
    except Exception as e:
        return f"<p style='color: red;'>Error visualizing SMILES: {e}</p>"



    
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def extract_protein_descriptors(sequence):
    try:
        analysis = ProteinAnalysis(sequence)
        descriptors = {
            "Protein Length": len(sequence),
            "Protein Molecular Weight": analysis.molecular_weight(),
            "Protein Aromaticity": analysis.aromaticity(),
            "Protein Instability Index": analysis.instability_index(),
            "Protein Isoelectric Point": analysis.isoelectric_point(),
            "Protein Gravy": analysis.gravy(),
            
        }

        # Amino acid composition
        amino_acid_percent = analysis.get_amino_acids_percent()
        for aa, percent in amino_acid_percent.items():
            descriptors[f"Protein Amino Acid Percent {aa}"] = percent

        return descriptors
    except Exception as e:
        print(f"Error processing protein sequence: {sequence}\n{e}")
        return {}




def PLInteraction(request):
    input1 = ''
    input2 = ''
    result = None
    smiles_visualization = None
    smiles_3D_visualization = None 
    ligand_descriptors= None
    protein_descriptors= None
    active_content_id = 'content-2'

    if request.method == 'POST':
        input1 = request.POST.get('input1', '')
        input2 = request.POST.get('input2', '')
        active_content_id = request.POST.get('activeContentId', 'content-2')
        result = f"Processed: SMILES = {input1}, Protein = {input2}"

        if input1:
            smiles_visualization = visualize_smiles(input1)
            smiles_3D_visualization=visualize_3D_smiles(input1)
            ligand_descriptors=generate_ligand_descriptors(input1)
            

        if input2:
            protein_descriptors=extract_protein_descriptors(input2)    

    return render(request, 'pl_int.html', {
        'input1': input1,
        'input2': input2,
        'result': result,
        'smiles_visualization': smiles_visualization,
        'smiles_3D_visualization':smiles_3D_visualization,
        'descriptors':ligand_descriptors,
        'active_content_id': active_content_id,
        'protein_descriptors':protein_descriptors,
    })


from django.core.files.storage import FileSystemStorage
from django.shortcuts import render





def process_input(request):
    if request.method == 'POST' and 'pdb_file' in request.FILES:
        # Handle the uploaded file
        pdb_file = request.FILES['pdb_file']
        fs = FileSystemStorage()  # Define where files are stored
        filename = fs.save(pdb_file.name, pdb_file)  # Save the file
        file_url = fs.url(filename)  # Get the file's URL if needed
        
        return render(request, 'pl_int.html', {
            'result_1': True,  # Pass this to indicate success
            'file_url': file_url,  # Optional: URL of the uploaded file
            'active_content_id': 'content-5',  # Preserve tab state
        })
    
    return render(request, 'pl_int.html', {
        'active_content_id': 'content-5',  # Preserve tab state
    })


'''
def sidebar_pl(request):
    
    if "username" not in request.session:
        messages.error(request, "Please log in first!")
        return redirect(f"{reverse('login')}?next={request.path}")
    
    username = request.session.get("username", None)  # Get username from session
    
    
    input1 = ''  # SMILES input (string)
    input2 = ''  # Protein input (possibly file path or string)
    input3=''
    input4=''
    db_name = request.POST.get("db", "default")
    result = None
    smiles_visualization = None
    smiles_3D_visualization = None
    ligand_descriptors = None
    protein_descriptors = None
    pdb_visual=None

    viewer_style="height:600px; width:700px; position:absolute; top:225px; left:1150px;"
    viewer_class = "viewer_3Dmoljs"
    data_pdb = ""  # Set your PDB data if available
    data_backgroundcolor= "0xffffff"
    data_style= "ballstick"
    data_ui= "true"
    

    #active_content_id = 'content-2'
    #uploaded_file_url = None
    #protein_descriptors_pdb=None
    #file_uploaded = False 

    if request.method == 'POST':
        # Handling non-file inputs (SMILES, Protein strings)
        input1 = request.POST.get('input1_pl', '')
        input2 = request.POST.get('input2_pl', '')
        input3 = request.FILES.get('input3_pl')
        input4 = request.FILES.get('input4_pl')

        if input3 is None:
           print(" Not Uploaded file:", input3)
        # Connect to MongoDB
        client = MongoClient(settings.DATABASES["mongo_db"]["CLIENT"]["host"])
        db = client[settings.DATABASES["mongo_db"]["NAME"]]  # Get the database
            
        # Collection (like a SQL table)
        collection = db.user_inputs  
            
        # Insert data
        collection.insert_one({"input1": input1, "input2": input2})
        
    
        result = f" Ki(nm) -> 7.5 "

        # Process SMILES and Protein inputs (if provided)
        if input1:
            smiles_visualization = visualize_smiles(input1)
            smiles_3D_visualization = visualize_3D_smiles(input1)
            ligand_descriptors = generate_ligand_descriptors(input1)

            
            # Connect to MongoDB
            client = MongoClient(settings.DATABASES["mongo_db"]["CLIENT"]["host"])
            db = client[settings.DATABASES["mongo_db"]["NAME"]]  # Select the DB
            collection = db.smiles_images  # Collection for storing SMILES images

            collection.insert_one({
                "image_data": smiles_visualization,  # Store the image in binary format
                "image_name": 'FirstImage',  # Store the image filename
                "uploaded_at": datetime.datetime.now()  # Optional: Store upload timestamp
            })

            



        if input2:
            protein_descriptors = extract_protein_descriptors(input2)

            

        if input3:
            pdb_visual= visualize_pdb(input3)

        if input4:
            pdb_visual= visualize_pdb(input4)

            

    return render(request, 'index_sidebar_pl.html', {
        'input1': input1,
        'input2': input2,
        'result': result,
        'smiles_visualization': smiles_visualization,
        'smiles_3D_visualization': smiles_3D_visualization,
        'descriptors': ligand_descriptors,
        'protein_descriptors': protein_descriptors,
        'username':username,
        'pdb_visual':pdb_visual,
        "viewer_style": viewer_style,
        "viewer_class": viewer_class,
        "data_pdb": data_pdb,
        "data_backgroundcolor": data_backgroundcolor,
        "data_style": data_style,
        "data_ui": data_ui,

    })       
'''

def sidebar_dd(request):

    if "username" not in request.session:
        messages.error(request, "Please log in first!")
        return redirect(f"{reverse('login')}?next={request.path}")
    
    username = request.session.get("username", None)  # Get username from session
    input1 = ''
    input2= ''
    result = None
    smiles_visualization_1 = None
    smiles_3D_visualization_1 = None 
    smiles_visualization_2 = None
    smiles_3D_visualization_2 = None 
    ligand_descriptors_1= None
    ligand_descriptors_2= None

    if request.method == 'POST':
        input1 = request.POST.get('input1_pl', '')
        input2 = request.POST.get('input2_pl', '')
        result = f" Interaction: Low "

        if input1:
            smiles_visualization_1 = visualize_smiles(input1)
            smiles_3D_visualization_1=visualize_3D_smiles(input1)
            ligand_descriptors_1=generate_ligand_descriptors(input1)

        if input2:
            smiles_visualization_2 = visualize_smiles(input2)
            smiles_3D_visualization_2=visualize_3D_smiles(input2)
            ligand_descriptors_2=generate_ligand_descriptors(input2)    

    return render(request, 'index_sidebar_dd.html', {
        'input1': input1,
        'input2': input2,
        'result': result,
        'smiles_visualization_1': smiles_visualization_1,
        'smiles_3D_visualization_1':smiles_3D_visualization_1,
        'smiles_visualization_2': smiles_visualization_2,
        'smiles_3D_visualization_2':smiles_3D_visualization_2,
        'descriptors_1':ligand_descriptors_1,
        'descriptors_2':ligand_descriptors_2,
        'username':username,
    })

def sidebar_pl(request):
    
    if "username" not in request.session:
        messages.error(request, "Please log in first!")
        return redirect(f"{reverse('login')}?next={request.path}")
    
    username = request.session.get("username", None)  # Get username from session

    # Initialize inputs and results
    input1 = ''  # SMILES input (string)
    input2 = ''  # Protein input (possibly file path or string)
    input3 = None
    input4 = None
    db_name = request.POST.get("db", "default")
    result = None
    smiles_visualization = None
    smiles_3D_visualization = None
    ligand_descriptors = None
    protein_descriptors = None
    pdb_visual = "Load the PDB file.."
    pdb_read="Load the PDB file.."

    viewer_style = "height:600px; width:700px; position:absolute; top:225px; left:1150px;"
    viewer_class = "viewer_3Dmoljs"
    data_pdb = ""  # Set your PDB data if available
    data_backgroundcolor = "0xffffff"
    data_style = "ballstick"
    data_ui = "true"

    if request.method == 'POST':
        # Get the hidden field to distinguish the forms.
        form_type = request.POST.get('form_type', 'project')  # default to project if not provided

        # For both forms, get the common fields (if any)
        input1 = request.POST.get('input1_pl', '')
        input2 = request.POST.get('input2_pl', '')
        # Files for project form and pdb form might be in different inputs
        input3 = request.FILES.get('input3_pl')  # From project form (optional)
        input4 = request.FILES.get('input4_pl')  # From PDB-only form

        # Connect to MongoDB (common for both)
        client = MongoClient(settings.DATABASES["mongo_db"]["CLIENT"]["host"])
        db = client[settings.DATABASES["mongo_db"]["NAME"]]
        collection = db.user_inputs  
        collection.insert_one({"input1": input1, "input2": input2})
        
        # Default result message (you can change this as needed)
        result = "7.5"

        if form_type == 'project':
            # Process the project form
            # SMILES and protein descriptors
            if input1:
                smiles_visualization = visualize_smiles(input1)
                smiles_3D_visualization = visualize_3D_smiles(input1)
                ligand_descriptors = generate_ligand_descriptors(input1)

                # Save SMILES image to MongoDB
                collection = db.smiles_images  # Collection for storing SMILES images
                collection.insert_one({
                    "image_data": smiles_visualization,
                    "image_name": 'FirstImage',
                    "uploaded_at": datetime.datetime.now()
                })

            if input2:
                protein_descriptors = extract_protein_descriptors(input2)
                

            # Process PDB file from project form if provided
            if input3:
                pdb_visual = visualize_pdb(input3)
                input2 = str(get_chain_a_sequence(input3))
                protein_descriptors = extract_protein_descriptors(input2)
                pdb_read=parse_atom_records_with_units_html(input3)
            else:
                print("Project form: Not Uploaded file for input3_pl")

        elif form_type == 'pdb':
            # Process the PDB file form
            if input4:
                pdb_visual= visualize_pdb(input4)
            

    context = {
        'input1': input1,
        'input2': input2,
        'result': result,
        'smiles_visualization': smiles_visualization,
        'smiles_3D_visualization': smiles_3D_visualization,
        'descriptors': ligand_descriptors,
        'protein_descriptors': protein_descriptors,
        'username': username,
        'pdb_visual': pdb_visual,
        "viewer_style": viewer_style,
        "viewer_class": viewer_class,
        "data_pdb": data_pdb,
        "data_backgroundcolor": data_backgroundcolor,
        "data_style": data_style,
        "data_ui": data_ui,
        "pdb_read":pdb_read
    }
    
    return render(request, 'index_sidebar_pl.html', context)




def sidebar_pp(request):
    if "username" not in request.session:
        messages.error(request, "Please log in first!")
        return redirect(f"{reverse('login')}?next={request.path}")
    
    username = request.session.get("username", None)  # Get username from session
    input1 = ''
    input2 = ''
    input3 = None
    input4 = None
    pdb_visual_1 = "Load the PDB-I file.."
    pdb_read_1="Load the PDB-I file.."
    pdb_visual_2 = "Load the PDB-II file.."
    pdb_read_2="Load the PDB-II file.."

    result = None
    smiles_visualization = None
    smiles_3D_visualization = None 
    protein_descriptors_1= None
    protein_descriptors_2= None

    viewer_style = "height:600px; width:700px; position:absolute; top:225px; left:1150px;"
    viewer_class = "viewer_3Dmoljs"
    data_pdb = ""  # Set your PDB data if available
    data_backgroundcolor = "0xffffff"
    data_style = "ballstick"
    data_ui = "true"
    

    if request.method == 'POST':
        input1 = request.POST.get('input1_pp', '')
        input2 = request.POST.get('input2_pp', '')
        input3 = request.FILES.get('input3_pp')  
        input4 = request.FILES.get('input4_pp')  
        
        result = f" Interaction: High"

        if input1:
            
            protein_descriptors_1=extract_protein_descriptors(input1)

        if input2:
            protein_descriptors_2=extract_protein_descriptors(input2) 

        if input3:
            pdb_visual_1 = visualize_pdb_pp(input3)
            input1 = str(get_chain_a_sequence(input3))
            protein_descriptors_1 = extract_protein_descriptors(input1)
            pdb_read_1= parse_atom_records_with_units_html(input3)
        else:
            print("Project form: Not Uploaded file for input3_pl")
            

        if input4:
            pdb_visual_2 = visualize_pdb_pp(input4)
            input2 = str(get_chain_a_sequence(input4))
            protein_descriptors_2 = extract_protein_descriptors(input2)
            pdb_read_2=parse_atom_records_with_units_html(input4)
        else:
            print("Project form: Not Uploaded file for input4_pl")

             




    return render(request, 'index_sidebar_pp.html', {
        'input1': input1,
        'input2': input2,
        'result': result,
        'smiles_visualization': smiles_visualization,
        'smiles_3D_visualization':smiles_3D_visualization,
        'protein_descriptors_1':protein_descriptors_1,
        'protein_descriptors_2':protein_descriptors_2,
        'username':username,
        "pdb_read_1":pdb_read_1,
        "pdb_read_2":pdb_read_2,
        "pdb_visual_1":pdb_visual_1,
        "pdb_visual_2":pdb_visual_2,
        "viewer_style": viewer_style,
        "viewer_class": viewer_class,
        "data_pdb": data_pdb,
        "data_backgroundcolor": data_backgroundcolor,
        "data_style": data_style,
        "data_ui": data_ui,
    })


def dropdown_mol_str(request):
    return render(request, 'dropdown_visualize_molstr.html')

def dropdown_prot_str(request):
    return render(request, 'dropdown_visualize_protstr.html')

def dropdown_li_desc(request):
    return render(request, 'dropdown_descriptor_li.html')

def dropdown_pro_desc(request):
    return render(request, 'dropdown_descriptor_pro.html')


from django.shortcuts import render, redirect
from django.contrib import messages
from django.contrib.auth.hashers import make_password, check_password
from pymongo import MongoClient

# Connect to MongoDB
client = MongoClient(settings.DATABASES["mongo_db"]["CLIENT"]["host"])
db = client[settings.DATABASES["mongo_db"]["NAME"]]  # Select the DB
users_collection = db["users"] # Collection for storing SMILES images


# Signup View
def signup_view(request):
    if "username" in request.session:
        return redirect("home")  # Redirect if already logged in
    if request.method == "POST":
        username = request.POST.get("username")
        email = request.POST.get("email")
        password1 = request.POST.get("password1")
        password2 = request.POST.get("password2")

        # Check if passwords match
        if password1 != password2:
            messages.error(request, "Passwords do not match!")
            return redirect("signup")

        # Check if username or email already exists
        existing_user = users_collection.find_one({"$or": [{"username": username}, {"email": email}]})
        if existing_user:
            messages.error(request, "Username or email already taken!")
            return redirect("signup")

        # Hash the password before saving
        hashed_password = make_password(password1)

        # Insert user into MongoDB
        users_collection.insert_one({
            "username": username,
            "email": email,
            "password": hashed_password  # Store the hashed password
        })

        messages.success(request, "Account created successfully! Please log in.")
        return redirect("login")

    return render(request, "signup.html")

from django.urls import reverse

# Login View
def login_view(request):
    if "username" in request.session:
        return redirect("home")  # Redirect if already logged in

    next_url = request.GET.get("next")  # Get 'next' parameter from URL

    if request.method == "POST":
        username = request.POST.get("username")
        password = request.POST.get("password")

        user = users_collection.find_one({"username": username})

        if user and check_password(password, user["password"]):
            request.session["username"] = user["username"]
            messages.success(request, "Login successful!")

            # Redirect back to the intended page or home if not set
            return redirect(next_url) if next_url else redirect("home")

        else:
            messages.error(request, "Invalid username or password!")
            return redirect(reverse("login") + f"?next={next_url}" if next_url else reverse("login"))

    return render(request, "login.html", {"next": next_url})

# Logout View
def logout_view(request):
    if "username" in request.session:
        del request.session["username"]  # Clear session
        messages.success(request, "Logged out successfully!")
    return redirect("home")


def visualize_pdb(pdbfile):
    """
    # Read the PDB file content
    with open(pdbfile, "r") as pdb_file:
        pdb_data = pdb_file.read()
    """

    pdb_data = pdbfile.read().decode('utf-8')

    # Create a py3Dmol viewer
    view = py3Dmol.view(width=800, height=600)

    # Load the PDB structure
    view.addModel(pdb_data, "pdb")  

    # Apply a cartoon representation
    view.setStyle({"cartoon": {"color": "spectrum"}})

    # Zoom to fit the structure
    view.zoomTo()

    embed_code = view._make_html()  

    return embed_code


def visualize_pdb_pp(pdbfile):
    """
    # Read the PDB file content
    with open(pdbfile, "r") as pdb_file:
        pdb_data = pdb_file.read()
    """

    pdb_data = pdbfile.read().decode('utf-8')

    # Create a py3Dmol viewer
    view = py3Dmol.view(width=800, height=300)

    # Load the PDB structure
    view.addModel(pdb_data, "pdb")  

    # Apply a cartoon representation
    view.setStyle({"cartoon": {"color": "spectrum"}})

    # Zoom to fit the structure
    view.zoomTo()

    embed_code = view._make_html()  

    return embed_code



def get_chain_a_sequence(pdb_input):
    """
    Extracts the primary sequence for Chain A from a PDB file using SEQRES records.
    
    Args:
        pdb_input (str or file-like): Either a file path to the PDB file or a file-like object.
        
    Returns:
        str: The one-letter FASTA sequence for chain A. Returns an empty string if not found.
    """
    # Mapping from three-letter to one-letter amino acid codes
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    residues_chain_a = []
    
    # If pdb_input is a string (file path), open it; otherwise assume it's a file-like object.
    if isinstance(pdb_input, (str, bytes)):
        # Open using the file path
        with open(pdb_input, "r") as pdb_file:
            lines = pdb_file.readlines()
    else:
        # pdb_input is assumed to be a file-like object (e.g., InMemoryUploadedFile)
        pdb_input.seek(0)  # ensure we're at the beginning
        # Read and decode the content (assuming it's UTF-8 encoded)
        content = pdb_input.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8")
        lines = content.splitlines()
    
    for line in lines:
        if line.startswith("SEQRES"):
            # Chain identifier is in column 12 (index 11)
            chain_id = line[11].strip()
            if chain_id == 'A':
                # Residues start at column 19 (index 18) and are space-separated
                residues = line[18:].split()
                residues_chain_a.extend(residues)
    
    # Convert three-letter codes to one-letter codes
    chain_a_sequence = "".join(three_to_one.get(res.upper(), 'X') for res in residues_chain_a)
    return chain_a_sequence


from tabulate import tabulate
import os

def parse_atom_records_with_units_html(pdb_input):
    """
    Parses ATOM records from a PDB file and returns a structured HTML table,
    including the dimensions/units in the headers.
    
    Args:
        pdb_input (str or file-like object): A path to the PDB file or a file-like object.
    
    Returns:
        str: A formatted HTML table of ATOM records.
    """
    atom_data = []
    
    # Determine if pdb_input is a file path or a file-like object.
    if isinstance(pdb_input, (str, os.PathLike)):
        with open(pdb_input, "r") as pdb_file:
            lines = pdb_file.readlines()
    else:
        # Assume file-like object, ensure reading from the beginning.
        pdb_input.seek(0)
        content = pdb_input.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8")
        lines = content.splitlines()
    
    for line in lines:
        if line.startswith("ATOM"):
            # Using fixed-width slicing based on PDB format specifications.
            record_name  = line[0:6].strip()
            serial       = line[6:11].strip()
            atom_name    = line[12:16].strip()
            alt_loc      = line[16].strip()
            residue_name = line[17:20].strip()
            chain_id     = line[21].strip()
            residue_seq  = line[22:26].strip()
            insertion    = line[26].strip()
            x            = line[30:38].strip()  # Coordinates in Angstroms (Å)
            y            = line[38:46].strip()  # Coordinates in Angstroms (Å)
            z            = line[46:54].strip()  # Coordinates in Angstroms (Å)
            occupancy    = line[54:60].strip()
            temp_factor  = line[60:66].strip()  # Temperature factor in Å²
            element      = line[76:78].strip()
            charge       = line[78:80].strip()
            
            atom_data.append([
                record_name, serial, atom_name, residue_name,
                chain_id, residue_seq, x, y, z, occupancy, temp_factor,
                element, charge
            ])
    
    # Include units/dimensions in the headers.
    headers = [
        "Rec", "Serial", "Atom", "Res", "Chain", "ResSeq", 
        "X (Å)", "Y (Å)", "Z (Å)", "Occupancy", "TempFactor (Å²)",
        "Elem", "Charge"
    ]
    # Return the table formatted as HTML
    return tabulate(atom_data, headers=headers, tablefmt="html")



def pl_results(ligand, sequence):
    pass


def dd_results(ligand, sequence):
    pass


def pp_results(ligand, sequence):
    pass