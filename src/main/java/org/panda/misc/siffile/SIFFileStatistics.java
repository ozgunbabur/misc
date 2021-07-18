package org.panda.misc.siffile;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class SIFFileStatistics
{
	public static void main(String[] args) throws IOException
	{
		printStatesForFile("/Users/ozgun/Documents/Analyses/CPTAC-LSCC-3.2/tumor-vs-normal/causative.sif");
	}

	public static void printStatesForFile(String file) throws IOException
	{
		Set<String> sources = new HashSet<>();
		Set<String> targets = new HashSet<>();
		Set<String> allRels = new HashSet<>();
		Map<String, Set<String>> typeToRels = new HashMap<>();

		Files.lines(Paths.get(file)).map(l -> l.split("\t"))/*.filter(t -> t.length > 3 && !t[3].isEmpty() && !t[3].equals("null"))*/.forEach(t ->
		{
			if (t.length > 2)
			{
				sources.add(t[0]);
				targets.add(t[2]);

				if (!typeToRels.containsKey(t[1])) typeToRels.put(t[1], new HashSet<>());

				String rel = t[0] + "\t" + t[1] + "\t" + t[2];
				typeToRels.get(t[1]).add(rel);
				allRels.add(rel);
			}
		});

		System.out.println("sources = " + sources.size());
		System.out.println("targets = " + targets.size() + "\n");
		System.out.println("interactions = " + allRels.size() + "\n");
		for (String type : typeToRels.keySet())
		{
			System.out.println(type + " = " + typeToRels.get(type).size());
		}
	}
}
