package org.panda.misc;

import org.biopax.paxtools.PaxtoolsMain;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.XReferrable;
import org.biopax.paxtools.model.level3.Xref;
import org.biopax.paxtools.pattern.miner.*;
import org.biopax.paxtools.pattern.util.Blacklist;
import org.panda.resource.HGNC;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class SIFGenerator
{
	public static void main(String[] args) throws Exception
	{
		useSIFSearcher();
	}

	public static void useSIFSearcher() throws Exception
	{
		String dir = "/Users/ozgun/Downloads/";
		SimpleIOHandler io = new SimpleIOHandler(BioPAXLevel.L3);
		Model model = io.convertFromOWL(new FileInputStream(dir + "electron-transport.owl"));

//		BlacklistGenerator bg = new BlacklistGenerator();
//		Blacklist blacklist = bg.generateBlacklist(model);
//		blacklist.write(dir + "blacklist-Reactome.txt");

		SIFSearcher ss = new SIFSearcher(ele -> {
			Set<String> ids = new HashSet<>();
			if (ele instanceof XReferrable)
			{
				for (Xref xref : ((XReferrable) ele).getXref())
				{
					if (xref.getDb().equals("hgnc symbol"))
					{
						ids.add(xref.getId());
					}
				}
			}
			System.out.println("ids = " + ids);
			return ids;
		}, SIFEnum.CATALYSIS_PRECEDES);

//		ss.setBlacklist(blacklist);
		Set<SIFInteraction> sifs = ss.searchSIF(model);
		System.out.println("sifs.size() = " + sifs.size());
		BufferedWriter writer = FileUtil.newBufferedWriter(dir + "electron-transport-catalysis-precedes.sif");
		sifs.forEach(sif -> FileUtil.writeln(sif.toString(), writer));
		writer.close();
	}

	public static void usePaxtoolsMain() throws Exception
	{
		String dir = "/home/babur/Documents/PC/";
//		String dir = "/home/babur/Documents/Temp/";

		PaxtoolsMain.toSifnx(new String[]{"toSIF", dir +

			"PathwayCommons9.All.BIOPAX.owl", dir + "PC.sif",
//		    "temp.owl", dir + "temp.sif",

			"-extended", "seqDb=hgnc", "chemDb=chebi"});
	}
}
