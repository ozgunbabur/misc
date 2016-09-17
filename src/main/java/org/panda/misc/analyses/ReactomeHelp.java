package org.panda.misc.analyses;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.misc.MutexReader;
import org.panda.misc.MutexReader.*;
import org.panda.utility.graph.Graph;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class ReactomeHelp
{
	public static final String MUTEX_DIR = "/home/babur/Documents/DARPA/BigMech/mutex/LGG/REACH-PC2v8";

	public static void main(String[] args)
	{
		printDifferentialEdgesAndPMCIDs();
	}

	static void printDifferentialEdgesAndPMCIDs()
	{
		Graph reach = new Graph("REACH", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		reach.load("/home/babur/Documents/DARPA/BigMech/mutex/networks/REACH.sif",
			Collections.emptySet(), Collections.singleton(reach.getEdgeType()));

		Graph pc = new Graph("PC", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		pc.load("/home/babur/Documents/DARPA/BigMech/mutex/networks/PC2v8.sif",
			Collections.emptySet(), Collections.singleton(reach.getEdgeType()));

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
}
