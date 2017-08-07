package org.panda.misc;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

/**
 * Created by babur on 3/3/16.
 */
public class GraphComparison
{
	public static void main(String[] args)
	{
		DirectedGraph g1 = new DirectedGraph("PC", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		g1.load("/home/babur/Documents/DARPA/BigMech/REACH-mutex/network/PC2v8.sif", new HashSet<>
			(Arrays.asList(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag())));
		DirectedGraph g2 = new DirectedGraph("REACH-reduced", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		g2.load("/home/babur/Documents/DARPA/BigMech/REACH-mutex/network/REACH-reduced.sif", new HashSet<>
			(Arrays.asList(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag())));
		DirectedGraph g3 = new DirectedGraph("REACH", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		g3.load("/home/babur/Documents/DARPA/BigMech/REACH-mutex/network/REACH.sif", new HashSet<>
			(Arrays.asList(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag())));

		g1.printVennIntersections(g2, g3);
	}
}
