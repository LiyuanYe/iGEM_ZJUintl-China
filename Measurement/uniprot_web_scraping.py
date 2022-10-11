"""
MAKE SURE YOUR ENVIRONMENT HAS THE FOLLOWING PACKAGES INSTALLED
"""
from urllib.request import urlopen
from urllib.error import HTTPError
from urllib.error import URLError
from bs4 import BeautifulSoup
import selenium
from selenium import webdriver
# DRIVER INCLUDED IN PACKAGE
# OR DOWNLOAD FROM https://chromedriver.chromium.org/downloads AND CHANGE PATH BELOW
from webdriver_manager.chrome import ChromeDriverManager

driver = webdriver.Chrome(ChromeDriverManager().install())
# driver = webdriver.Chrome("D:\\Jetbrains\\PyCharmPro\\projects\\SRTP21\\chromedriver.exe")
import time


"""
FUNCTIONS
"""
# read potential markers file
def read_markers(filepath):
    with open(filepath, 'r') as infile:
        lines = infile.readlines()
        potential_markers = [line.split("\n")[0] for line in lines]
        infile.close()
    return potential_markers

# get page function
def get_page_by_url(url):
    try:
        # html = urlopen(url)
        # soup = BeautifulSoup(html.read(), 'html.parser')
        driver.get(url)
        content = driver.page_source
        soup = BeautifulSoup(content, 'html.parser')
        return soup
    except HTTPError as e:
        print(e)
    except URLError as e:
        print("Server not found!")
    else:
        print("Connection Successful...")

# get entry id function
def uniprot_search_gene_by_name_species(target_gene, target_species):
    page = get_page_by_url("https://www.uniprot.org/uniprot/?query=" + target_gene.lower())
    try:
        res_table = page.find_all('table', {'class': 'data-table'})[0]
    except IndexError:
        return -1
    tbody = res_table.tbody
    res_entries = []
    for tr in tbody.children:
        # get gene name
        # gene_name = tr.find_all('span', {'class': 'shortName'})
        all_cols = tr.find_all('span')
        gene_name = all_cols[4]
        if gene_name:
            gene_name = gene_name.get_text()
        # get species name
        species = tr.find_all('a')[1]
        if species:
            species = species.get_text()
        # get reviewed
        reviewed = tr.find_all('span', {'class': 'entry-title__status icon--reviewed'})
        if gene_name and species:
            name_species_match = target_species in species and target_gene.lower() in gene_name.lower()
        else:
            name_species_match = False
        if name_species_match and reviewed:
            res_entries.append(tr.find_all('a')[0].text)
        elif name_species_match and len(res_entries) == 0:
            print(gene_name + "only has unreviewed item(s)!")
    return res_entries

# check membrane protein function - DON'T USE THIS FUNCTION!
def check_membrane_protein_by_uniprot(entry):
    on_membrane = False
    driver.get("https://www.uniprot.org/uniprot/" + entry.lower())
    # expand normally invisible shadow root using javascript
    try:
        location = driver.execute_script("return document.querySelector('sib-swissbiopics-sl').shadowRoot.getElementById('swissbiopic')")
    except selenium.common.exceptions.JavascriptException:
        on_membrane = False
        return on_membrane
    # we can continue if present
    outer_html = location.get_attribute('outerHTML')
    location_soup = BeautifulSoup(outer_html, 'html.parser')
    swissBioPicsSlData = location_soup.find_all('div', {'id': 'swissBioPicsSlData'})
    if swissBioPicsSlData:
        for isoform in swissBioPicsSlData[0].children:
            has_plasma_membrane_item = isoform.find_all('li', {'class': 'Plasma_Membrane'})
            def has_other_location_membrane_item(isoform):
                # check
                secreted = isoform.find_all('li', {'class': 'Extracellular_region_or_secreted'})
                endo = isoform.find_all('li', {'class': 'Endoplasmic_Reticulum'})
                other_locations = isoform.find_all('a', {'class': 'subcell_name'})
                for location in other_locations:
                    if location.text == 'Membrane ' and not secreted and not endo:
                        return True
                return False
            if has_plasma_membrane_item or has_other_location_membrane_item(isoform):
                on_membrane = True
                break
    # if no UniProt term, we judge based on GO term
    else:
        go_annotation = location_soup.find_all('div', {'id': 'table-go_annotation'})
        for isoform in go_annotation[0].children:
            has_plasma_membrane_item = isoform.find_all('li', {'class': 'Plasma_Membrane'})
            if has_plasma_membrane_item:
                on_membrane = True
                break
    return on_membrane

# check membrane protein by new version uniprot keyword function - USE THIS ONE INSTEAD
def check_membrane_protein_by_uniprot_keyword(entry):
    on_membrane = False
    driver.get("https://www.uniprot.org/uniprotkb/" + entry.lower() + "/entry")
    time.sleep(3)
    content = driver.page_source
    soup = BeautifulSoup(content, 'html.parser')
    target_section = soup.find_all('section', {'id': 'subcellular_location'})
    if target_section:
        keywords = target_section[0].find_all('ul', {'class': 'expandable-list no-bullet'})
        if keywords:
            keywords = keywords[0]
        else:
            return on_membrane
    else:
        print("can't find subcell section")
    locations = [keyword.get_text() for keyword in keywords.children]
    if " #Cell membrane" in locations:
        on_membrane = True
    return on_membrane

# largest test if membrane related function
def uniprot_test_membrane_gene_by_name_species(target_gene, target_species):
    is_membrane_protein = False
    entries = uniprot_search_gene_by_name_species(target_gene, target_species)
    if entries == -1:
        return False
    for entry in entries:
        time.sleep(0.1)
        if check_membrane_protein_by_uniprot_keyword(entry):
            is_membrane_protein = True
    return is_membrane_protein

"""
TEST AREA
"""
inpath = 'D:\\SRTP2021\\web_scraping_test\\test.txt'
outpath = 'D:\\SRTP2021\\web_scraping_test\\out.txt'
markers = read_markers(inpath)

res_entries = uniprot_search_gene_by_name_species('cyba', 'Mus musculus')

print(check_membrane_protein_by_uniprot_keyword('Q61462')) # cyba

print(uniprot_test_membrane_gene_by_name_species('cdkn2a', 'Mus musculus'))

"""
MAIN PROCESS
"""
# inpath = 'D:\\SRTP2021\\potential_markers\\human_proteome_rna1.csv'
# outpath = 'D:\\SRTP2021\\potential_markers\\human_proteome_rna1_res.csv'

def main_function(inpath=None, outpath=None, object=None, mode='file', developer_mode=False, membrane_proteins_check=None):
    try:
        if mode == 'file':
            potential_markers = read_markers(inpath)
        elif mode == 'object':
            potential_markers = object
        # main programm, uses a list of genes
        cnt = 0
        cnt_membrane  = 0
        membrane_proteins = []
        for gene in potential_markers:
            print("Searching for ", gene, " in UniProt Database...")
            time.sleep(0.1)
            cnt += 1
            if uniprot_test_membrane_gene_by_name_species(gene, 'Mus musculus'):
                cnt_membrane += 1
                membrane_proteins.append(gene)
                print(gene, " is a membrane protein!")
                # developer mode
                if developer_mode:
                    if gene not in membrane_proteins_check:
                        input("This gene is not in membrane proteins, press ENTER to continue")
            else:
                # developer mode
                if developer_mode:
                    if gene in membrane_proteins_check:
                        input("This gene should be a membrane protein, press ENTER to continue")
                print(gene, " is not a membrane protein.")
            print(cnt, " of ", len(potential_markers), " checked... ", cnt_membrane," membrane proteins found\n")
        print("Process finished, found ", len(membrane_proteins), " membrane proteins.")
        # output result
        if mode == 'file':
            with open(outpath, 'a') as outfile:
                for protein in membrane_proteins:
                    outfile.write(protein + '\n')
                outfile.close()
    except Exception:
        print()
    finally:
        return membrane_proteins


test_membrane_found = main_function(inpath=inpath, outpath=outpath)

# membrane_proteins_found = main_function(inpath=inpath, outpath=outpath)


