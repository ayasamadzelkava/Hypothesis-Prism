import javafx.animation.AnimationTimer
import javafx.application.Application
import javafx.geometry.Point3D
import javafx.scene.*
import javafx.scene.input.MouseEvent
import javafx.scene.input.ScrollEvent
import javafx.scene.paint.Color
import javafx.scene.paint.PhongMaterial
import javafx.scene.shape.Cylinder
import javafx.scene.shape.Sphere
import javafx.scene.text.Text
import javafx.scene.transform.Rotate
import javafx.stage.Stage
import kotlin.math.acos

// -------------------------------------
// Tree Data Structures, Parser, and Layout
// -------------------------------------

data class TreeNode(
    val label: String,
    var level: Int,
    val type: String,
    var children: MutableList<TreeNode> = mutableListOf(),
    var position: Point3D = Point3D(0.0, 0.0, 0.0),
    var u: Double = 0.0  // Horizontal parameter for layout (0 = left, 1 = right)
)

fun parseTreeDefinition(defStr: String): TreeNode {
    val trimmedDef = if (defStr.contains(":")) defStr.substringAfter(":").trim() else defStr.trim()
    if (trimmedDef.isEmpty()) throw IllegalArgumentException("Tree definition string is empty.")
    val segments = trimmedDef.split(";").map { it.trim() }.filter { it.isNotEmpty() }
    if (segments.isEmpty()) throw IllegalArgumentException("No valid segments found.")
    val nodesMap = mutableMapOf<String, TreeNode>()
    val childrenMapping = mutableMapOf<String, List<String>>()
    for (segment in segments) {
        val parts = segment.split("->")
        val leftPart = parts[0].trim() // e.g., "A(lvl0, source)"
        val label = leftPart.substringBefore("(").trim()
        val inside = leftPart.substringAfter("(").substringBefore(")").trim()
        val tokens = inside.split(",").map { it.trim() }
        val level = tokens[0].substringAfter("lvl").toIntOrNull() ?: 0
        val type = if (tokens.size > 1) tokens[1] else ""
        nodesMap[label] = TreeNode(label, level, type)
        if (parts.size > 1) {
            val childrenPart = parts[1].trim()
            if (childrenPart.isNotEmpty()) {
                val childLabels = childrenPart.split(",").map { it.trim() }.filter { it.isNotEmpty() }
                childrenMapping[label] = childLabels
            }
        }
    }
    // Link children.
    for ((parentLabel, childLabels) in childrenMapping) {
        val parentNode = nodesMap[parentLabel]!!
        for (childLabel in childLabels) {
            val childNode = nodesMap.getOrPut(childLabel) { TreeNode(childLabel, 0, "") }
            parentNode.children.add(childNode)
        }
    }
    val rootLabel = segments[0].substringBefore("(").trim()
    return nodesMap[rootLabel]!!
}

fun assignLeafOrderAndPropagate(root: TreeNode) {
    val leaves = mutableListOf<TreeNode>()
    fun dfs(node: TreeNode) {
        if (node.children.isEmpty()) leaves.add(node) else node.children.forEach { dfs(it) }
    }
    dfs(root)
    val n = leaves.size
    leaves.forEachIndexed { index, node ->
        node.u = if (n > 1) index.toDouble() / (n - 1) else 0.5
    }
    fun propagate(node: TreeNode) {
        if (node.children.isNotEmpty()) {
            node.children.forEach { propagate(it) }
            node.u = node.children.map { it.u }.average()
        }
    }
    propagate(root)
}

fun computeMaxLevel(root: TreeNode): Int {
    var maxLevel = root.level
    for (child in root.children) {
        maxLevel = maxOf(maxLevel, computeMaxLevel(child))
    }
    return maxLevel
}

fun interpolate(p1: Point3D, p2: Point3D, t: Double): Point3D {
    return Point3D(
        p1.x + (p2.x - p1.x) * t,
        p1.y + (p2.y - p1.y) * t,
        p1.z + (p2.z - p1.z) * t
    )
}

fun layoutTree(root: TreeNode, source: Point3D, leftLeaf: Point3D, rightLeaf: Point3D) {
    val maxLevel = computeMaxLevel(root).toDouble()
    fun layoutNode(node: TreeNode) {
        val fraction = if (maxLevel == 0.0) 0.0 else node.level / maxLevel
        val baseline = interpolate(leftLeaf, rightLeaf, node.u)
        node.position = interpolate(source, baseline, fraction)
        node.children.forEach { layoutNode(it) }
    }
    layoutNode(root)
}

// -------------------------------------
// JavaFX Helpers for 3D Objects
// -------------------------------------

fun create3DLine(start: Point3D, end: Point3D, radius: Double, color: Color): Cylinder {
    val diff = end.subtract(start)
    val height = diff.magnitude()
    val mid = start.midpoint(end)
    val cylinder = Cylinder(radius, height)
    cylinder.material = PhongMaterial(color)
    val yAxis = Point3D(0.0, 1.0, 0.0)
    val diffNorm = diff.normalize()
    val angle = Math.toDegrees(acos(yAxis.dotProduct(diffNorm)))
    val axis = yAxis.crossProduct(diffNorm)
    if (axis.magnitude() > 0) cylinder.transforms.add(Rotate(angle, axis))
    cylinder.translateX = mid.x
    cylinder.translateY = mid.y
    cylinder.translateZ = mid.z
    return cylinder
}

fun Point3D.midpoint(other: Point3D): Point3D {
    return Point3D((this.x + other.x) / 2, (this.y + other.y) / 2, (this.z + other.z) / 2)
}

fun createVertexMarker(position: Point3D, radius: Double, color: Color): Sphere {
    val sphere = Sphere(radius)
    sphere.material = PhongMaterial(color)
    sphere.translateX = position.x
    sphere.translateY = position.y
    sphere.translateZ = position.z
    return sphere
}

fun createTextLabel(textStr: String, position: Point3D, color: Color): Text {
    val text = Text(textStr)
    text.fill = color
    text.translateX = position.x
    text.translateY = position.y
    text.translateZ = position.z
    text.scaleX = 2.0
    text.scaleY = 2.0
    return text
}

// -------------------------------------
// Graph Data Model and Builder (Combined IDgraph tree)
// -------------------------------------

data class GraphNode(
    val id: Int,
    val label: String,
    val position: Point3D,
    val treeType: String, // e.g., "I" or "D"
    val level: Int
)

data class GraphEdge(
    val from: Int,
    val to: Int,
    val color: Color = Color.PURPLE
)

data class Graph(
    val nodes: List<GraphNode>,
    val adjMatrix: Array<BooleanArray>,
    val edges: List<GraphEdge>
)

fun mapTreeNodesToGraphIds(
    node: TreeNode,
    treeType: String,
    mapping: MutableMap<TreeNode, Int>,
    nodeList: MutableList<GraphNode>
) {
    val id = nodeList.size
    mapping[node] = id
    nodeList.add(GraphNode(id, node.label, node.position, treeType, node.level))
    node.children.forEach { mapTreeNodesToGraphIds(it, treeType, mapping, nodeList) }
}

fun collectTreeEdges(
    node: TreeNode,
    mapping: Map<TreeNode, Int>,
    edgeList: MutableList<GraphEdge>,
    edgeColor: Color
) {
    for (child in node.children) {
        edgeList.add(GraphEdge(mapping[node]!!, mapping[child]!!, edgeColor))
        collectTreeEdges(child, mapping, edgeList, edgeColor)
    }
}

fun collectLeaves(node: TreeNode, leaves: MutableList<TreeNode>) {
    if (node.children.isEmpty()) leaves.add(node)
    else node.children.forEach { collectLeaves(it, leaves) }
}

fun buildIDGraph(iTreeRoot: TreeNode, dTreeRoot: TreeNode): Graph {
    val graphNodes = mutableListOf<GraphNode>()
    val graphEdges = mutableListOf<GraphEdge>()
    val mappingI = mutableMapOf<TreeNode, Int>()
    val mappingD = mutableMapOf<TreeNode, Int>()

    // Map nodes from I tree (we tag them as "I") and from D tree (tagged as "D")
    mapTreeNodesToGraphIds(iTreeRoot, "I", mappingI, graphNodes)
    mapTreeNodesToGraphIds(dTreeRoot, "D", mappingD, graphNodes)

    // Add intra-tree parent–child edges.
    collectTreeEdges(iTreeRoot, mappingI, graphEdges, Color.GREEN)
    collectTreeEdges(dTreeRoot, mappingD, graphEdges, Color.BLUE)

    // Collect leaves from both trees.
    val iLeaves = mutableListOf<TreeNode>()
    val dLeaves = mutableListOf<TreeNode>()
    collectLeaves(iTreeRoot, iLeaves)
    collectLeaves(dTreeRoot, dLeaves)

    // Add all–to–all connection edges between I leaves and D leaves.
    for (iLeaf in iLeaves) {
        for (dLeaf in dLeaves) {
            val from = mappingI[iLeaf]!!
            val to = mappingD[dLeaf]!!
            graphEdges.add(GraphEdge(from, to, Color.PURPLE))
        }
    }

    // Build an adjacency matrix.
    val n = graphNodes.size
    val matrix = Array(n) { BooleanArray(n) { false } }
    for (edge in graphEdges) {
        matrix[edge.from][edge.to] = true
        matrix[edge.to][edge.from] = true
    }

    return Graph(graphNodes, matrix, graphEdges)
}

fun buildGraphGroup(graph: Graph, nodeRadius: Double): Group {
    val group = Group()
    // Create markers and labels.
    for (node in graph.nodes) {
        // Use different colors for I and D nodes.
        val color = if (node.treeType == "I") Color.GREEN else Color.BLUE
        val marker = createVertexMarker(node.position, nodeRadius, color)
        val label = createTextLabel(node.label, node.position, Color.BLACK)
        group.children.addAll(marker, label)
    }
    // Draw edges.
    for (edge in graph.edges) {
        val nodeFrom = graph.nodes[edge.from]
        val nodeTo = graph.nodes[edge.to]
        val line = create3DLine(nodeFrom.position, nodeTo.position, 0.5, edge.color)
        line.viewOrder = -1.0
        group.children.add(line)
    }
    return group
}

// -------------------------------------
// Main Application: Combined IDgraph Tree
// -------------------------------------

class GraphTetrahedronApp : Application() {

    override fun start(primaryStage: Stage) {
        // Define tetrahedron vertices.
        // We use these points to position the trees on different faces.
        val iSource = Point3D(0.0, 100.0, 0.0)            // I tree source
        val dSource = Point3D(-100.0, -100.0, 100.0)         // D tree source
        val dLeftEndpoint = Point3D(100.0, -100.0, 100.0)     // D tree left endpoint
        val dRightEndpoint = Point3D(0.0, -100.0, -100.0)     // D tree right endpoint

        // For the I tree, compute endpoints for its leaf baseline.
        val iLeftEndpoint = Point3D(
            iSource.x + 0.5 * (dLeftEndpoint.x - iSource.x),
            iSource.y + 0.5 * (dLeftEndpoint.y - iSource.y),
            iSource.z + 0.5 * (dLeftEndpoint.z - iSource.z)
        )
        val iRightEndpoint = Point3D(
            iSource.x + 0.5 * (dRightEndpoint.x - iSource.x),
            iSource.y + 0.5 * (dRightEndpoint.y - iSource.y),
            iSource.z + 0.5 * (dRightEndpoint.z - iSource.z)
        )

        // Build tetrahedron edges and vertex markers.
        val tetraEdgeColor = Color.LIGHTGRAY.deriveColor(0.0, 1.0, 1.0, 0.3)
        val edge1 = create3DLine(iSource, dSource, 1.0, tetraEdgeColor)
        val edge2 = create3DLine(iSource, dLeftEndpoint, 1.0, tetraEdgeColor)
        val edge3 = create3DLine(iSource, dRightEndpoint, 1.0, tetraEdgeColor)
        val edge4 = create3DLine(dSource, dLeftEndpoint, 1.0, tetraEdgeColor)
        val edge5 = create3DLine(dSource, dRightEndpoint, 1.0, tetraEdgeColor)
        val edge6 = create3DLine(dLeftEndpoint, dRightEndpoint, 1.0, tetraEdgeColor)
        val edgesGroup = Group(edge1, edge2, edge3, edge4, edge5, edge6)
        edgesGroup.viewOrder = 0.0

        val vertexMarkers = Group(
            createVertexMarker(iSource, 5.0, Color.GREEN),
            createVertexMarker(dSource, 5.0, Color.ORANGE),
            createVertexMarker(dLeftEndpoint, 5.0, Color.ORANGE),
            createVertexMarker(dRightEndpoint, 5.0, Color.ORANGE)
        )
        val vertexLabels = Group(
            createTextLabel("I", iSource, Color.GREEN),
            createTextLabel("d", dSource, Color.ORANGE),
            createTextLabel("m1", dLeftEndpoint, Color.ORANGE),
            createTextLabel("m2", dRightEndpoint, Color.ORANGE)
        )

        // Define tree specification strings.
        val iTreeDefinition = "A(lvl0, source)->B,C; B(lvl1, int)->D,E; C(lvl1, int)->F; D(lvl2, int)->; E(lvl2,int)->; F(lvl2,int)->"
        val dTreeDefinition = "X(lvl0, source)->Y,Z; Y(lvl1, int)->P; Z(lvl1, int)->Q,R; Q(lvl2,measurement)->; R(lvl2,measurement)->; P(lvl2,measurement)->"

        // Parse the trees.
        val iTreeRoot = parseTreeDefinition(iTreeDefinition)
        val dTreeRoot = parseTreeDefinition(dTreeDefinition)
        assignLeafOrderAndPropagate(iTreeRoot)
        assignLeafOrderAndPropagate(dTreeRoot)
        layoutTree(iTreeRoot, iSource, iLeftEndpoint, iRightEndpoint)
        layoutTree(dTreeRoot, dSource, dLeftEndpoint, dRightEndpoint)

        // Build the combined graph.
        val idGraph = buildIDGraph(iTreeRoot, dTreeRoot)
        val idGraphGroup = buildGraphGroup(idGraph, 3.0)

        // Combine tetrahedron edges, vertex markers, labels, and the IDgraph.
        val sceneGroup = Group(edgesGroup, vertexMarkers, vertexLabels, idGraphGroup)

        // Center the tetrahedron.
        val center = Point3D(
            (iSource.x + dSource.x + dLeftEndpoint.x + dRightEndpoint.x) / 4,
            (iSource.y + dSource.y + dLeftEndpoint.y + dRightEndpoint.y) / 4,
            (iSource.z + dSource.z + dLeftEndpoint.z + dRightEndpoint.z) / 4
        )
        sceneGroup.translateX = -center.x
        sceneGroup.translateY = -center.y
        sceneGroup.translateZ = -center.z

        // Apply rotations.
        val rotateX = Rotate(20.0, Rotate.X_AXIS)
        val rotateY = Rotate(-20.0, Rotate.Y_AXIS)
        sceneGroup.transforms.addAll(rotateX, rotateY)

        // Set up camera and scene.
        val camera = PerspectiveCamera(true)
        camera.translateZ = -800.0
        camera.nearClip = 0.1
        camera.farClip = 2000.0

        val rootGroup = Group(sceneGroup)
        val scene = Scene(rootGroup, 800.0, 600.0, true, SceneAntialiasing.BALANCED)
        scene.fill = Color.WHITE
        scene.camera = camera

        // Interactivity: Mouse drag to rotate, scroll to zoom.
        var anchorX = 0.0
        var anchorY = 0.0
        var anchorAngleX = 0.0
        var anchorAngleY = 0.0

        scene.addEventHandler(MouseEvent.MOUSE_PRESSED) { event ->
            anchorX = event.sceneX
            anchorY = event.sceneY
            anchorAngleX = rotateX.angle
            anchorAngleY = rotateY.angle
        }
        scene.addEventHandler(MouseEvent.MOUSE_DRAGGED) { event ->
            rotateX.angle = anchorAngleX - (event.sceneY - anchorY)
            rotateY.angle = anchorAngleY + (event.sceneX - anchorX)
        }
        scene.addEventHandler(ScrollEvent.SCROLL) { event ->
            camera.translateZ += event.deltaY
        }

        // Animation timer to keep text labels facing the camera.
        object : AnimationTimer() {
            override fun handle(now: Long) {
                sceneGroup.children.filterIsInstance<Text>().forEach { text ->
                    text.transforms.setAll(
                        Rotate(-rotateY.angle, Rotate.Y_AXIS),
                        Rotate(-rotateX.angle, Rotate.X_AXIS)
                    )
                }
            }
        }.start()

        primaryStage.title = "Combined IDgraph Tree"
        primaryStage.scene = scene
        primaryStage.show()
    }
}

fun main(args: Array<String>) {
    Application.launch(GraphTetrahedronApp::class.java, *args)
}
