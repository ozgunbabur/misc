package org.panda.misc;

import org.panda.utility.CollectionUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Anil uses a different format for RPPA platform files. Because that it is easier to process rows in R when rows
 * contain only one gene symbol, they duplicate rows when
 *
 * @author Ozgun Babur
 */
public class MDAndersonPlatformToCausalPathPlatformFormat
{
	public static void main(String[] args) throws IOException
	{
//		checkOldFileMatching();
//		convert("/home/babur/Documents/Analyses/SMMART/Patient1/RPPA/antibodyMap.txt",
//			"/home/babur/Documents/Analyses/SMMART/Patient1/RPPA/mdanderson-platform.txt");
	}

	public static void convert(String mdFile, String cpFile) throws IOException
	{
		Map<String, Map<String, String>> siteMap = new HashMap<>();
		Map<String, String> effectMap = new HashMap<>();

		Files.lines(Paths.get(mdFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[2];
			String sym = t[3];
			String sites = t[4];
			String effect = t[5];

			if (!siteMap.containsKey(id)) siteMap.put(id, new HashMap<>());
			Map<String, String> map = siteMap.get(id);
			if (!map.containsKey(id))
			{
				map.put(sym, sites);
				effectMap.put(id, effect);
			}
			else
			{
				if (map.get(sym).length() < sites.length())
				{
					map.put(sym, sites);
				}

				if (!effectMap.get(id).equals(effect))
				{
					System.err.println("Conflicting effects: " + id);
				}
			}
		});

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(cpFile));
		writer.write("ID\tSymbols\tSites\tEffect");

		for (String id : siteMap.keySet())
		{
			Map<String, String> map = siteMap.get(id);

			List<String> symList = map.keySet().stream().sorted().collect(Collectors.toList());

			String aSite = map.values().iterator().next();
			if (aSite.equals("NA") || aSite.isEmpty())
			{
				writer.write("\n" + id + "\t" + CollectionUtil.merge(symList, " ") + "\t\t");
			}
			else
			{
				List<String> siteList = map.keySet().stream().sorted().map(map::get).collect(Collectors.toList());
				writer.write("\n" + id + "\t" + CollectionUtil.merge(symList, " ") + "\t" +
					CollectionUtil.merge(siteList, " ") + "\t" + effectMap.get(id));
			}
		}

		writer.close();
	}

	public static void checkOldFileMatching() throws IOException
	{
		Map<String, String> oldMap = new HashMap<>();
		Files.lines(Paths.get("/home/babur/Documents/Analyses/JQ1/platform.txt")).skip(1).forEach(l ->
		{
			String id = l.split("\t")[1].replaceAll("p", "_p").toUpperCase().replaceAll("_P", "_p");
			oldMap.put(id, id + l.substring(l.indexOf("\t", l.indexOf("\t") + 1)));
		});

		Map<String, String> newMap = new HashMap<>();
		Files.lines(Paths.get("/home/babur/Documents/Analyses/SMMART/Patient1/RPPA/mdanderson-platform.txt")).skip(1).forEach(l ->
		{
			String id = l.split("\t")[0].replaceAll("p", "_p").toUpperCase().replaceAll("_P", "_p");
			newMap.put(id, id + l.substring(l.indexOf("\t")));
		});

		String[] header = Files.lines(Paths.get("/home/babur/Documents/Analyses/SMMART/Patient1/RPPA/TCGA-BRCA-L4.csv"))
			.findFirst().get().split(",");

		BufferedWriter writer = Files.newBufferedWriter(Paths.get("/home/babur/Documents/Analyses/SMMART/Patient1/RPPA/tcga-platform.txt"));
		writer.write("ID\tSymbols\tSites\tEffect");

		for (int i = 4; i < header.length; i++)
		{
			if (oldMap.containsKey(header[i]))
			{
				writer.write("\n" + oldMap.get(header[i]));
			}
			else if (newMap.containsKey(header[i]))
			{
				writer.write("\n" + newMap.get(header[i]) + "\tattention");
			}
			else
			{
				writer.write("\n" + header[i] + "\t\t\t\tfill");
			}
		}

		writer.close();
	}

	public static void convertMDAndersonValuesFile()
	{

	}
}
