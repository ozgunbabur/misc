package org.panda.misc;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.network.IPTMNet;
import org.panda.resource.network.PathwayCommons;
import org.panda.resource.network.PhosphoNetworks;
import org.panda.resource.network.SignedPC;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;

import java.util.*;

/**
 * Created by babur on 3/3/16.
 */
public class GraphComparison
{
	public static void main(String[] args)
	{

	}

	public static void compareSignedDirected()
	{
		// load PC
		Map<SignedType, DirectedGraph> pcGraphs = SignedPC.get().getAllGraphs();
		Set<String> pc = new HashSet<>();
		pcGraphs.values().forEach(graph -> pc.addAll(getRelsAsSet(graph)));

		// load PhosphoNetworks
		Set<String> phosNet = getRelsAsSet(PhosphoNetworks.get().getGraph());

		// load IPTMNet
		Set<String> iptmNet = getRelsAsSet(IPTMNet.get().getGraph());

	}

	private static Set<String> getRelsAsSet(DirectedGraph graph)
	{
		Set<String> set = new HashSet<>();
		graph.getOneSideSymbols(true).stream().forEach(s ->
			graph.getDownstream(s).forEach(t -> set.add(s + "\t" + graph.getEdgeType() + "\t" + t)));
		return set;
	}

	public static void compareBasicSIF()
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
