package org.panda.misc;

import org.biopax.paxtools.PaxtoolsMain;

/**
 * @author Ozgun Babur
 */
public class SIFGenerator
{
	public static void main(String[] args) throws Exception
	{
		String dir = "/home/babur/Documents/PC/";
//		String dir = "/home/babur/Documents/Temp/";

		PaxtoolsMain.toSifnx(new String[]{"toSIF", dir +

			"PathwayCommons9.All.BIOPAX.owl", dir + "PC.sif",
//		    "temp.owl", dir + "temp.sif",

			"-extended", "seqDb=hgnc", "chemDb=chebi"});
	}
}
