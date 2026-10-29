🧬 Bioinformatics 3D Project
Explore protein structures in 3D, inspect protein sequences, and experiment with a structure-prediction workflow.
![🚀 Open Live Demo](https://img.shields.io/badge/🚀_Open-Live_Demo-FF4B4B?style=for-the-badge)
![Open in Colab](https://img.shields.io/badge/Google_Colab-Open_Prediction_Notebook-F9AB00?style=for-the-badge&logo=googlecolab&logoColor=white)
![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat&logo=python&logoColor=white)
![Streamlit](https://img.shields.io/badge/Built_with-Streamlit-FF4B4B?style=flat&logo=streamlit&logoColor=white)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
> **Try it online:** [Launch the Bioinformatics 3D Project](https://bioinformatics-3d-project-ai.streamlit.app/)  
> **Structure-prediction workflow:** [Open the Google Colab notebook](https://colab.research.google.com/github/sagnikdas1610/bioinformatics-3d-project/blob/main/AI_Structure_Generator.ipynb)
---
📌 Overview
The Bioinformatics 3D Project is an educational web application for exploring protein structures and examining protein sequence properties. It combines an interactive Streamlit interface with a separate Google Colab notebook for experimenting with protein-structure prediction using ESMFold.
The project is designed for students, learners, and anyone interested in structural bioinformatics. Upload a supported protein structure, explore its 3D representation, and inspect sequence-level measurements in one workspace.
✨ Features
Interactive 3D protein viewer — inspect uploaded protein structures in the browser.
PDB structure upload — load supported `.pdb` files into the viewer.
Protein sequence input — enter an amino-acid sequence in raw or FASTA format.
Sequence analysis — examine supported protein properties and amino-acid composition.
Structure exploration controls — use the molecular representations and coloring options available in the app.
Prediction workflow — open the companion Google Colab notebook to run the separate structure-prediction workflow.
Browser-based interface — access the deployed Streamlit app without installing the project locally.
Feature availability depends on the current deployed version and the input structure.
🖼️ Screenshots
Screenshots help visitors understand the project before opening it. Add screenshots captured from the actual running application; do not use mockups as if they were real application output.
Create a folder named `docs/screenshots/` in this repository, then upload these images:
Screenshot	Suggested filename	What to capture
Main dashboard	`dashboard.png`	Homepage and overall layout
3D protein viewer	`protein-viewer.png`	A loaded structure, such as hemoglobin
Sequence analysis	`sequence-analysis.png`	Sequence metrics and composition results
Upload workflow	`upload-workflow.png`	PDB upload and controls
After uploading the images, uncomment the image lines below (remove the `<!--` and `-->` markers):
<!--
### Main dashboard
![Bioinformatics 3D Project dashboard](docs/screenshots/dashboard.png)

### Interactive protein viewer
![Interactive 3D protein viewer](docs/screenshots/protein-viewer.png)

### Protein sequence analysis
![Protein sequence analysis](docs/screenshots/sequence-analysis.png)

### Structure upload workflow
![PDB upload workflow](docs/screenshots/upload-workflow.png)
🧰 Technology Stack
Technology	Purpose
Python	Application logic and analysis
Streamlit	Interactive web interface and deployment
Biopython	Protein sequence analysis utilities
py3Dmol	Browser-based molecular visualization
Google Colab / ESMFold workflow	Separate protein-structure prediction experiment
🚀 Use the Live App
Open the live website.
Upload a supported `.pdb` structure file.
If you want to analyze sequence properties, paste the matching amino-acid sequence into the sequence field.
Explore the 3D structure and the analysis results provided by the app.
Example structure files
Download experimental structures from the RCSB Protein Data Bank:
1CRN — Crambin
1UBQ — Ubiquitin
1LYZ — Lysozyme
4HHB — Hemoglobin
Use the sequence belonging to the chain you intend to analyze. A multi-chain PDB file may contain multiple protein chains; their sequences should not be combined unless the analysis specifically requires it.
🧪 Run Locally
Requirements
Python installed on your computer
Git (optional, for cloning the repository)
An internet connection to install dependencies
Installation
```bash
git clone https://github.com/sagnikdas1610/bioinformatics-3d-project.git
cd bioinformatics-3d-project
python -m pip install -r requirements.txt
```
Start the Streamlit app
```bash
python -m streamlit run app.py
```
Streamlit will print a local URL, usually `http://localhost:8501`. Open that address in your browser.
🧠 Protein-Structure Prediction Workflow
The repository includes `AI_Structure_Generator.ipynb`, which can be opened in Google Colab.
![Open in Google Colab](https://colab.research.google.com/assets/colab-badge.svg)
Prediction is a separate workflow from the deployed visualization and sequence-analysis website. Successful prediction depends on the notebook, available compute resources, model/service availability, and the sequence supplied.
⚠️ Scope and Limitations
This is an educational project, not a substitute for specialist structural-bioinformatics validation.
Visualizing an experimental PDB structure is different from predicting a new structure.
A predicted structure should not be treated as experimentally confirmed.
Results depend on the input sequence, PDB contents, and current app implementation.
Large structures may require more time and browser resources to render.
🔭 Future Improvements
Add more sample structures and a clearer example gallery.
Improve multi-chain selection and chain-specific sequence extraction.
Add downloadable analysis reports and richer structure metadata.
Provide clearer prediction status, confidence information, and error handling.
Add automated tests and documented example outputs.
Expand mobile responsiveness and accessibility.
🤝 Contributions and Feedback
Suggestions, bug reports, and improvements are welcome. Open an issue in this repository with a clear description of the problem or proposed enhancement.
👨‍💻 Author
Sagnik Das  
Computer Science & Engineering student interested in Python, AI, and bioinformatics.
GitHub: @sagnikdas1610
Project repository: bioinformatics-3d-project
Live demo: Open the web app
📄 License
This project is released under the MIT License.
