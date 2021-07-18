package org.panda.misc;

import org.panda.resource.MSigDB;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class ExamineSIFCoverageOfGeneSets
{
	public static void main(String[] args)
	{
		Map<String, Set<String>> geneSets = MSigDB.get().getSetsNameFiltered(name -> name.startsWith("BIOCARTA_") || name.startsWith("KEGG_") || name.startsWith("PID_") || name.startsWith("REACTOME_") || name.startsWith("WP_"));
		List<String[]> rels = loadSIF("/Users/ozgun/Documents/Analyses/CausalPath-data/causal-priors.txt");
		printCoverage(rels, geneSets);
	}

	private static List<String[]> loadSIF(String filename)
	{
		return FileUtil.lines(filename).map(l -> l.split("\t")).filter(t -> t.length > 2)
			.map(t -> new String[]{t[0], t[2]}).collect(Collectors.toList());
	}

	private static void printCoverage(List<String[]> rels, Map<String, Set<String>> geneSets)
	{
		Set<String> relGenes = rels.stream().map(s -> s[0]).collect(Collectors.toSet());
		relGenes.addAll(rels.stream().map(s -> s[1]).collect(Collectors.toSet()));

		double totalGeneSets = geneSets.size();

		long hit1 = geneSets.keySet().stream().filter(k -> !CollectionUtil.getIntersection(geneSets.get(k), relGenes).isEmpty()).count();

		long hit2 = geneSets.keySet().stream().filter(k -> rels.stream().anyMatch(r -> geneSets.get(k).contains(r[0]) && geneSets.get(k).contains(r[1]))).count();

		System.out.println("totalGeneSets = " + totalGeneSets);
		System.out.println("Gene coverage     = " + (hit1 / totalGeneSets));
		System.out.println("Relation coverage = " + (hit2 / totalGeneSets));
	}
}
