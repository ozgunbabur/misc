package org.panda.misc.analyses;

import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.Named;
import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.misc.MutexReader;
import org.panda.misc.MutexReader.*;
import org.panda.resource.Pubmed;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class ReactomeHelp
{
	public static final String MUTEX_DIR = "/home/babur/Documents/DARPA/BigMech/mutex/LGG/REACH-PC2v8";

	public static void main(String[] args) throws IOException
	{
//		printDifferentialEdgesAndPMCIDs();
//		printRhoGEFRelations();
//		printRhoGEFRelationsVersion2();
//		printOverlap();

		printPhosphoRelsTextEvidence();
	}

	static void printDifferentialEdgesAndPMCIDs()
	{
		DirectedGraph reach = new DirectedGraph("REACH", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		reach.load("/home/babur/Documents/DARPA/BigMech/mutex/networks/REACH.sif",
			Collections.singleton(reach.getEdgeType()));

		DirectedGraph pc = new DirectedGraph("PC", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		pc.load("/home/babur/Documents/DARPA/BigMech/mutex/networks/PC2v8.sif",
			Collections.singleton(reach.getEdgeType()));

		Set<String> remember = new HashSet<>();

		Set<Group> groups = MutexReader.readMutexResults(MUTEX_DIR);
		groups.stream().filter(g -> g.score < 0.05).forEach(g ->
			g.genes.forEach(g1 -> g.genes.stream().filter(g2 -> !g2.equals(g1)).forEach(g2 ->
			{
				String key = g1 + "\t" + reach.getEdgeType() + "\t" + g2;
				if (!remember.contains(key))
				{
					remember.add(key);
					Set<String> rDwn = reach.getDownstream(g1);
					Set<String> pDwn = pc.getDownstream(g1);

					if (rDwn.contains(g2) && !pDwn.contains(g2))
					{
						System.out.println(key + "\t" + reach.getMediatorsInString(g1, g2));
					}
				}
			})));
	}

	static void printRhoGEFRelations() throws IOException
	{
		String dir = "/home/babur/Documents/Analyses/Reactome-help/";
		Set<String> reacSet = new HashSet<>();

		Set<String> rhos = new HashSet<>();
		Set<String> gefs = new HashSet<>();
		Files.lines(Paths.get(dir + "RHO_GTPase_GEFs_Aug10_2017.txt")).skip(1)
			.map(l -> l.split("\t")).forEach(t ->
		{
			rhos.add(t[0]);
			gefs.add(t[1]);
			if (t.length > 4 && !t[4].isEmpty())
			{
				reacSet.add(t[1] + " " + t[0]);
			}
		});

		Set<String> reachSet = new HashSet<>();
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + "GEF-Rho-in-Reach.txt"));
		writer.write("Source\tRelation-type\tTarget\tProcess-IDs\tLiterature-reference");
		Files.lines(Paths.get("/home/babur/Documents/DARPA/BigMech/REACH-mutex/network/REACH.sif")).forEach(l ->
		{
			String[] t = l.split("\t");
			if (t[1].startsWith("controls") && gefs.contains(t[0]) && rhos.contains(t[2]))
			{
				FileUtil.lnwrite(l, writer);
				reachSet.add(t[0] + " " + t[2]);
			}
		});
		writer.close();

		CollectionUtil.printNameMapping("Reactome", "Reach");
		CollectionUtil.printVennCounts(reacSet, reachSet);
	}

	static void printRhoGEFRelationsVersion2() throws IOException
	{
		Map<String, String> reach = new HashMap<>();
		Files.lines(Paths.get("/home/babur/Documents/DARPA/BigMech/REACH-wo-generics.sif")).map(l -> l.split("\t"))
			.forEach(t ->
			{
				String key = t[0] + " " + t[2];
				reach.put(key, t[4]);
			});

		String dir = "/home/babur/Documents/Analyses/Reactome-help/";
		String fileToFix = dir + "GEF-Rho-in-Reach_PMID_Reach_Reactome_Comparison-cleaned.txt";

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + "GEF-Rho_Reach_Reactome_Comparison.txt"));
//		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(System.out));
		writer.write(Files.lines(Paths.get(fileToFix)).findFirst().get());

		Files.lines(Paths.get(fileToFix)).skip(1)
			.map(l -> l.split("\t")).forEach(t ->
		{
			String key = t[1] + " " + t[0];
			FileUtil.lnwrite(t[0] + "\t" + t[1] + "\t" + t[2] + "\t", writer);
			System.out.print("\n" + t[0] + "\t" + t[1] + "\t" + t[2] + "\t");

			if (reach.containsKey(key))
			{
				String s = reach.get(key);
				String pmids = extractPMIDS(s);
				String text = replacePMCIDS(s);
				FileUtil.write(pmids + "\t" + text, writer);
				System.out.print(pmids + "\t" + text);
			}
			else
			{
				FileUtil.write("\t", writer);
				System.out.print("\t");
			}

			FileUtil.write("\t" + (t.length > 5 ? t[5] : ""), writer);
			System.out.print("\t" + (t.length > 5 ? t[5] : ""));
		});

	}

	static String extractPMIDS(String s)
	{
		List<String> list = new ArrayList<>();
		String[] t = s.split(";");

		for (String tt : t)
		{
			String pmc = tt.substring(0, tt.indexOf(":"));
			assert pmc.startsWith("PMC");
			String pm = Pubmed.get().getPMID(pmc);
			if (pm != null) list.add("PMID" + pm);
		}

		return CollectionUtil.merge(list, ";");
	}

	static String replacePMCIDS(String s)
	{
		List<String> list = new ArrayList<>();
		String[] t = s.split(";");

		for (String tt : t)
		{
			String pmc = tt.substring(0, tt.indexOf(":"));
			assert pmc.startsWith("PMC");
			String pm = Pubmed.get().getPMID(pmc);
			if (pm != null) list.add("PMID" + pm +  tt.substring(tt.indexOf(":")));
		}

		return CollectionUtil.merge(list, ";");
	}

	static void printOverlap() throws IOException
	{
		Set<String> reach = new HashSet<>();
		Set<String> react = new HashSet<>();
		int[] c = new int[]{0};

		Files.lines(Paths.get("/home/babur/Documents/Analyses/Reactome-help/GEF-Rho_Reach_Reactome_Comparison.txt"))
			.map(l -> l.split("\t")).forEach(t ->
		{
			String key = t[0] + t[1];
			if (t.length > 2 && !t[2].isEmpty() && !t[2].startsWith("N/D")) react.add(key);
			if (t.length > 3 && !t[3].isEmpty())
			{
				reach.add(key);

				if (t[2].startsWith("NO") && !t[2].contains(";")) c[0]++;
			}
		});

		CollectionUtil.printNameMapping("Reactome", "Reach");
		CollectionUtil.printVennCounts(react, reach);
		System.out.println("\n\n" + c[0]);
	}

	static void printPhosphoRelsTextEvidence() throws IOException
	{
		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer/BigMech/Diff-view/";
		SimpleIOHandler io = new SimpleIOHandler();
		Model model = io.convertFromOWL(new FileInputStream(dir + "model.owl"));

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + "reach-and-pc-with-text.txt"));

		Files.lines(Paths.get(dir + "reach-and-pc.sif")).map(l -> l.split("\t")).forEach(t ->
		{
			Set<String> locs = Arrays.stream(t[4].split(";")).map(s -> s.substring(1)).collect(Collectors.toSet());

			Set<String> set = new HashSet<>();
			for (String id : t[3].split(" "))
			{
				if (id.startsWith("Control") || id.startsWith("Catalysis"))
				{
					Named ele = (Named) model.getByID(id);

					for (String com : ele.getComment())
					{
						if (locs.stream().filter(com::contains).findAny().isPresent())
						{
							set.add(com);
						}
					}
				}
			}

			FileUtil.writeln(ArrayUtil.getString("\t", t[0], t[1], t[2], t[4], CollectionUtil.merge(set, ";")), writer);
		});

		writer.close();
	}
}
