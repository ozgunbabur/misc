package org.panda.misc.bigmech;

import org.panda.misc.MutexReader;
import org.panda.utility.graph.DirectedGraph;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * This class is for printing the relations that REACH detects, but not PC, and those relations are in mutex relation
 * in a cancer study, detected when REACH network was used as the network resource.
 *
 * @author Ozgun Babur
 */
public class REACHDiffForMutex
{
	public static final String BASE = "/media/babur/6TB1/REACH-mutex/";
	public static final String NETWORK_DIR = BASE + "network/";
	public static final String REACH_NETWORK_FILE = NETWORK_DIR + "REACH.sif";
	public static final String PC_NETWORK_FILE = NETWORK_DIR + "PC2v8.sif";
	public static final String RUN_DIR = BASE + "run/";
	public static final String REACH_DIR = "/REACH-PC2v8";
	public static final String PC_DIR = "/PC2v8";
	public static final String[] STUDY = new String[]{"BRCA", "GBM", "LGG"};
	public static String REL_TYPE = "controls-state-change-of";
	public static final double SCORE_THR = 0.05;

	public static void main(String[] args) throws IOException
	{
		for (String study : STUDY)
		{
			printForStudy(study);
		}
	}

	static void printForStudy(String study) throws IOException
	{
		System.out.println("study = " + study);
		Map<String, Set<String>> reach = readMutexNetwork(RUN_DIR + study + REACH_DIR, SCORE_THR);
		Map<String, Set<String>> pc = readMutexNetwork(RUN_DIR + study + PC_DIR, SCORE_THR);

		removeMap2FromMap1(reach, pc);

		DirectedGraph rNetwork = readGraph(REACH_NETWORK_FILE);
		DirectedGraph pNetwork = readGraph(PC_NETWORK_FILE);
		Map<String, Map<String, String>> det = readREACHDetails();

		List<String> result = new ArrayList<>();

		for (String g1 : reach.keySet())
		{
			for (String g2 : reach.get(g1))
			{
				if (rNetwork.getDownstream(g1).contains(g2) && !pNetwork.getDownstream(g1).contains(g2))
				{
					result.add(det.get(g1).get(g2));
				}
				if (rNetwork.getDownstream(g2).contains(g1) && !pNetwork.getDownstream(g2).contains(g1))
				{
					result.add(det.get(g2).get(g1));
				}
			}
		}
		Collections.sort(result);
		result.forEach(System.out::println);
	}

	static Map< String, Set<String>> readMutexNetwork(String dir, double thr)
	{
		Map< String, Set<String>> map = new HashMap<>();
		MutexReader.readMutexResults(dir).stream().filter(g -> g.getScore() <= thr).forEach(g ->
		{
			for (String g1 : g.genes)
			{
				for (String g2 : g.genes)
				{
					if (g1.compareTo(g2) < 0)
					{
						if (!map.containsKey(g1)) map.put(g1, new HashSet<>());
						map.get(g1).add(g2);
					}
				}
			}
		});

		return map;
	}

	static void removeMap2FromMap1(Map< String, Set<String>> map1, Map< String, Set<String>> map2)
	{
		Set<String> rem = new HashSet<>();
		for (String gene : map1.keySet())
		{
			if (map2.containsKey(gene)) map1.get(gene).removeAll(map2.get(gene));
			if (map1.get(gene).isEmpty()) rem.add(gene);
		}
		for (String gene : rem)
		{
			map1.remove(gene);
		}
	}

	static Map<String, Map<String, String>> readREACHDetails() throws IOException
	{
		Map<String, Map<String, String>> map = new HashMap<>();

		Files.lines(Paths.get(REACH_NETWORK_FILE)).forEach(l ->
		{
			String[] token = l.split("\t");
			if (token[1].equals(REL_TYPE))
			{
				if (!map.containsKey(token[0])) map.put(token[0], new HashMap<>());
				map.get(token[0]).put(token[2], l);
			}
		});

		return map;
	}

	static DirectedGraph readGraph(String file)
	{
		DirectedGraph graph = new DirectedGraph();
		graph.load(file, Collections.singleton(REL_TYPE));
		return graph;
	}
}
