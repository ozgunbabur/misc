package org.panda.misc.analyses;

import org.panda.resource.HGNC;
import org.panda.resource.MGI;
import org.panda.resource.SiteMappingMouseToHuman;
import org.panda.resource.network.UniProt;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class Satpathy
{
	public static final String DIR = "/Users/ozgun/Documents/Analyses/Satpathy/";

	public static void main(String[] args) throws IOException
	{
		prepareDataFiles();
	}

	public static void prepareDataFiles() throws IOException
	{
		String outFile = DIR + "data.txt";
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("ID\tSymbols\tSites\tEffect\tM/L\tH/L");

		Files.lines(Paths.get(DIR + "proteome.csv")).skip(23).map(l -> l.split("\t")).forEach(t ->
		{
			Set<String> mSet = new HashSet<>(Arrays.asList(t[2].split(";")));
			List<String> hList = getHumanSymbols(mSet);

			if (!hList.isEmpty())
			{
				FileUtil.lnwrite(CollectionUtil.merge(hList, "_") + "\t" + CollectionUtil.merge(hList, " ") + "\t\t", writer);

				FileUtil.tab_write(t[9], writer);
				FileUtil.tab_write(t[11], writer);
			}
		});

		SiteMappingMouseToHuman mth = new SiteMappingMouseToHuman();


		Files.lines(Paths.get(DIR + "phosphoproteome.csv")).skip(30).map(l -> l.split("\t")).forEach(t ->
		{
			Set<String> mSet = new HashSet<>(Arrays.asList(t[4].split(";")));
			Set<String> hSet = new HashSet<>(getHumanSymbols(mSet));

			if (hSet.isEmpty()) return;

			List<String> mUPs = Arrays.asList(t[0].split(";"));
			String aa = t[11];

			assert aa.length() == 1;

			List<String> mSites = Arrays.asList(t[1].substring(1).split(";"));
			Map<String, String> hMap = new HashMap<>();

			for (int i = 0; i < mUPs.size(); i++)
			{
				if (!mUPs.get(i).contains("-"))
				{
					String mUP = mUPs.get(i);
					String mSite = aa + mSites.get(i);

					Map<String, List<String>> converted = mth.mapToHumanSite(mUP, mSite);

					Map<String, String> convSymMap = new HashMap<>();

					for (String hUP : converted.keySet())
					{
						String hSym = HGNC.get().getSymbol(hUP);
						if (hSym != null) convSymMap.put(hSym, converted.get(hUP).get(0));
					}

					for (String hSym : convSymMap.keySet())
					{
						if (hSet.contains(hSym) && !hMap.containsKey(hSym))
						{
							hMap.put(hSym, convSymMap.get(hSym));
						}
					}
				}
			}

			if (!hMap.isEmpty())
			{
				double log2 = Math.log(2);

				String id = "";
				String syms = "";
				String sites = "";

				for (String sym : hMap.keySet())
				{
					String site = hMap.get(sym);

					id += sym + "_" + site + "_";
					syms += sym + " ";
					sites += site + " ";
				}

				id = id.substring(0, id.length() - 1);
				syms = syms.trim();
				sites = sites.trim();

				FileUtil.lnwrite(id + "\t" + syms + "\t" + sites + "\t", writer);
				FileUtil.tab_write(Math.log(Double.valueOf(t[17].isEmpty()?"1":t[17])) / log2, writer);
				FileUtil.tab_write(Math.log(Double.valueOf(t[18].isEmpty()?"1":t[18])) / log2, writer);
			}
		});


		writer.close();
	}

	private static List<String> getHumanSymbols(Collection<String> mSyms)
	{
		Set<String> hSet = new HashSet<>();
		for (String mSym : mSyms)
		{
			Set<String> hums = MGI.get().getCorrespondingHumanSymbols(mSym);
			hSet.addAll(hums);
		}
		List<String> hSyms = new ArrayList<>(hSet);
		hSyms.sort(String::compareTo);
		return hSyms;
	}
}
