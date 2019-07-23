package org.panda.misc.analyses;

import org.panda.resource.HGNC;
import org.panda.resource.MGI;
import org.panda.resource.SiteMappingMouseToHuman;
import org.panda.resource.UniProtSequence;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author Ozgun Babur
 */
public class PaulBPhosphoproteomics
{
	public static final String DIR = "/home/babur/Documents/Analyses/PaulBarnes/";

	public static void main(String[] args) throws IOException
	{
		convertMouseToHuman();
	}

	public static void convertMouseToHuman() throws IOException
	{
		int tnum = 3;
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "T" + tnum + "/data.txt"));
		writer.write("ID\tGenes\tSites\tEffect\tFold-change");

		Files.lines(Paths.get(DIR + "Treatment" + tnum + ".csv")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String pep = t[0];
			String markedPep = t[1];
			String name = t[2];
			double fc = Double.valueOf(t[3]);

			int start = getStartPosition(pep, name);
			if (start < 1) System.err.println("Cannot find peptide: " + name);

			List<String> sites = getModifications(markedPep, start);
			String upID = UniProtSequence.get().getIDOfName(name);

			Map<String, String> hSiteMap = new HashMap<>();

			Map<String, List<String>> mapping = SiteMappingMouseToHuman.get().mapToHumanSite(upID, sites.toArray(new String[sites.size()]));

			String mSym = MGI.get().getSymbol(upID);
			Set<String> hSyms = MGI.get().getCorrespondingHumanSymbols(mSym);

			for (String hUP : mapping.keySet())
			{
				String sym = HGNC.get().getSymbol(hUP);
				if (sym != null && hSyms.contains(sym))
				{
					hSiteMap.put(sym, CollectionUtil.merge(mapping.get(hUP), "|"));
				}
			}
			if (!hSiteMap.isEmpty())
			{
				// todo continue

				String syms = CollectionUtil.merge(hSiteMap.keySet().stream().sorted().collect(Collectors.toList()), " ");
				String siteStr = CollectionUtil.merge(hSiteMap.keySet().stream().sorted().map(hSiteMap::get).collect(Collectors.toList()), " ");

				String id = CollectionUtil.merge(hSiteMap.keySet().stream().sorted().map(sym -> sym + "-" + hSiteMap.get(sym).replaceAll("\\|", "-")).collect(Collectors.toList()), "-");

				FileUtil.lnwrite(id + "\t" + syms + "\t" + siteStr + "\t\t" + fc,
					writer);
			}
		});

		writer.close();
	}

	static int getStartPosition(String peptide, String uniprotName)
	{
		return UniProtSequence.get().getStartLocation(uniprotName, peptide);
	}

	static List<String> getModifications(String markedPeptide, int startPos)
	{
		markedPeptide = markedPeptide.substring(2, markedPeptide.length() - 2);
		List<String> sites = new ArrayList<>();
		String[] t = markedPeptide.split("\\*");
		for (int i = 0; i < t.length - 1; i++)
		{
			char lett = t[i].charAt(t[i].length() - 1);
			int pos = startPos + t[i].length() - 1;

			sites.add(lett + "" + pos);

			startPos += t[i].length();
		}

		return sites;
	}


}
