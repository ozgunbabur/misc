package org.babur.misc;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.Graph;

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
		Graph g1 = new Graph("PC", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		g1.load("/home/babur/Documents/mutex/networks/PC2v7.sif", Collections.<String>emptySet(), new HashSet<>
			(Arrays.asList(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag(), SIFEnum.CONTROLS_EXPRESSION_OF.getTag())));
		Graph g2 = new Graph("Fries290K", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		g2.load("/home/babur/Documents/mutex/networks/fries290K.sif", Collections.<String>emptySet(), new HashSet<>
			(Arrays.asList(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag(), SIFEnum.CONTROLS_EXPRESSION_OF.getTag())));
		Graph g3 = new Graph("Leidos", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		g3.load("/home/babur/Documents/mutex/networks/leidos.sif", Collections.<String>emptySet(), new HashSet<>
			(Arrays.asList(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag(), SIFEnum.CONTROLS_EXPRESSION_OF.getTag())));
		Graph g4 = new Graph("fries-Ras", "is-upstream-of");
		g4.load("/home/babur/Documents/mutex/networks/fries1K.sif", Collections.<String>emptySet(), new HashSet<>
			(Arrays.asList("is-upstream-of")));

		g1.printVennIntersections(g2, g3, g4);
//		g2.printVennIntersections(g4);
	}
}
