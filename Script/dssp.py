import sys
def extraction_atome(fichier_pdb, fichier_text):
    """ Lit un fichier pdb et renvoie les coordonées de tous les atomes
        impliqués dans la liaison peptidique (C,O,H,N).

    Parametres
    ----------
    fichier_pdb: fichier pdb contenant les coordonnées des atomes
    fichier_text: fichier de sortie qui contient que les coordonnées
                des atomes impliqués dans la liaison hydrogene

    Returns
    -------
    position: liste du numero des differents residus
    fichier_text

    """
    modele = 1
    atome = 0
    nom = ["C","O","N","H"]
    with open(fichier_pdb, "r") as pdb, open(fichier_text, "w") as text:
        for ligne in pdb:
            if ligne.startswith("ATOM") and ligne[12:16].strip() in nom:
                atome += 1
                posi=ligne[24:27].strip()
                nam=ligne[12:16].strip()
                resid=ligne[17:20].strip()
                x = ligne[32:38]
                y = ligne[40:46]
                z = ligne[48:54]
                text.write("{}	{}	{}	{}	{}	{} \n".format(posi,resid,nam,x, y, z))
    return fichier_text

