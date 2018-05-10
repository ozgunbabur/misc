package org.panda.misc;

import org.panda.resource.PCPathwayHGNC;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class SIFPathwayGroupsAdder
{
	public static void main(String[] args) throws IOException
	{
		String sifWoExt = "/home/babur/Documents/RPPA/TCGA/PNNL/correlation-based-phospho/causative";
		addGroups(sifWoExt + ".sif", sifWoExt + "-with-pathways.sif", 0.1);
		Runtime.getRuntime().exec("cp " + sifWoExt + ".format " + sifWoExt + "-with-pathways.format");
	}

	public static void addGroups(String inSIF, String outSIF, double pvalThr) throws IOException
	{
		List<String> lines = Files.lines(Paths.get(inSIF)).collect(Collectors.toList());

		Set<String> genes = lines.stream().map(l -> l.split("\t"))
			.map(t -> t.length == 1 ? new String[]{t[0]} : new String[]{t[0], t[2]})
			.flatMap(Arrays::stream).distinct().collect(Collectors.toSet());

		List<Group> groups = new ArrayList<>();

		while (true)
		{
			PCPathwayHGNC pcp = new PCPathwayHGNC();
			Map<String, Double> pvals = pcp.calculateEnrichment(genes, 3, 250);
			String groupID = getMin(pvals, pvalThr, pcp, 3, genes);
			if (groupID == null) break;

			Set<PCPathwayHGNC.Pathway> pathways = pcp.getPathways(groupID);
			Group group = new Group(pathways, genes);
			groups.add(group);

			genes.removeAll(group.genes);
			if (genes.size() < 2) return;
		}

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outSIF));
		lines.forEach(l -> FileUtil.writeln(l, writer));
		groups.forEach(g -> FileUtil.writeln(g.toString(), writer));
		writer.close();
	}


	private static String getMin(Map<String, Double> map, double pvalThr, PCPathwayHGNC pcp, int minMember,
		Set<String> queryGenes)
	{
		String s = null;
		double min = 1;

		for (String id : map.keySet())
		{
			Double val = map.get(id);
			if (val < min && val <= pvalThr)
			{
				Set<PCPathwayHGNC.Pathway> pathways = pcp.getPathways(id);
				int intersect = CollectionUtil.countOverlap(pathways.iterator().next().getGenes(), queryGenes);

				if (intersect >= minMember)
				{
					s = id;
					min = val;
				}
			}
		}

		return s;
	}

	static class Group
	{
		String name;
		String tooltip;
		Set<String> genes;

		public Group(Set<PCPathwayHGNC.Pathway> pathways, Set<String> freeGenes)
		{
			name = getName(pathways);
			tooltip = getTooltip(pathways);
			genes = new HashSet<>(pathways.iterator().next().getGenes());
			genes.retainAll(freeGenes);
		}

		private static String getName(Set<PCPathwayHGNC.Pathway> pathways)
		{
			String shortest = null;
			int size = Integer.MAX_VALUE;

			for (PCPathwayHGNC.Pathway pathway : pathways)
			{
				if (pathway.getName().length() < size)
				{
					shortest = pathway.getName();
					size = shortest.length();
				}
			}
			return shortest;
		}

		private static String getTooltip(Set<PCPathwayHGNC.Pathway> pathways)
		{
			StringBuilder sb = new StringBuilder();

			for (PCPathwayHGNC.Pathway pathway : pathways)
			{
				sb.append(pathway.getSource()).append(": ").append(pathway.getName()).append("\n");
			}

			return sb.toString().trim();
		}

		@Override
		public String toString()
		{
			return "$group$\t" + name + "\t" + CollectionUtil.merge(genes, " ") + "\t" + tooltip;
		}
	}
}
