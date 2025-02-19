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
// Unified Collapse Rendering on IDgraph (Unified Graph) Based on Internal Levels
// -------------------------------------
/*
  This function performs collapse on the unified IDgraph tree.
  Each node has an "internal level" computed as follows:
    - For nodes from the I tree: internal level = declared level.
    - For nodes from the D tree: internal level = (n + m) – declared level,
      where n = (# levels in I tree) and m = (# levels in D tree).

  The function traverses the unified tree structure (which is built by merging the two trees
  along their intra-tree edges) and if a node’s internal level is in the collapse set, it is hidden.
  In addition, if after processing its children no visible descendant is found, then the node itself
  is treated as a "bridge" (i.e. returned as a visible leaf) so that connecting (cross) edges can be drawn.

  Note: This function is applied separately to the I and D subtrees (with appropriate internal level mapping)
  and then the cross edges are drawn between the visible leaves of each. This produces a unified IDgraph collapse.
*/
fun renderCollapsedIDGraph(
    root: TreeNode,
    effectiveParent: TreeNode?,
    collapseLevels: Set<Int>,
    getInternalLevel: (TreeNode) -> Int,
    nodeColor: Color,
    nodeRadius: Double
): Pair<Group, List<TreeNode>> {
    val group = Group()
    val visibleLeaves = mutableListOf<TreeNode>()
    val internalLevel = getInternalLevel(root)
    val isCollapsed = internalLevel in collapseLevels

    // Draw node if not collapsed.
    if (!isCollapsed) {
        group.children.add(createVertexMarker(root.position, nodeRadius, nodeColor))
        group.children.add(createTextLabel(root.label, root.position, Color.BLACK))
    }
    // In unified collapse, if the node is collapsed, we do NOT use its position as an anchor.
    val newEffective = if (!isCollapsed) root else effectiveParent

    // Process children.
    for (child in root.children) {
        val childInternal = getInternalLevel(child)
        val childCollapsed = childInternal in collapseLevels
        // If a child is visible, draw an edge from newEffective to child.
        if (newEffective != null && !childCollapsed) {
            group.children.add(create3DLine(newEffective.position, child.position, 1.0, nodeColor))
        }
        val (childGroup, childLeaves) = renderCollapsedIDGraph(child, newEffective, collapseLevels, getInternalLevel, nodeColor, nodeRadius)
        group.children.add(childGroup)
        visibleLeaves.addAll(childLeaves)
    }
    // If no visible leaves were found in the subtree and this node is not collapsed,
    // treat this node as a visible leaf (bridge) so that cross edges can attach.
    if (visibleLeaves.isEmpty() && !isCollapsed) {
        visibleLeaves.add(root)
    }
    return Pair(group, visibleLeaves)
}

// -------------------------------------
// Main Application: Unified IDgraph Tree with Collapse on All Levels
// -------------------------------------

class GraphTetrahedronApp : Application() {

    // Global collapse set: any internal level in this set will be collapsed.
    private val globalCollapseLevels: MutableSet<Int> = mutableSetOf()

    // Group holding the rendered (collapsed) unified IDgraph tree.
    private var idGraphCollapseGroup: Group = Group()

    // Single overlay group for toggle buttons.
    private var toggleOverlayGroup: Group = Group()

    // Global references to the parsed I and D trees.
    private lateinit var iTreeRoot: TreeNode
    private lateinit var dTreeRoot: TreeNode

    // The overall scene group.
    private lateinit var sceneGroup: Group

    // For mouse rotation.
    private var anchorX = 0.0
    private var anchorY = 0.0
    private var anchorAngleX = 0.0
    private var anchorAngleY = 0.0
    private val rotateX = Rotate(20.0, Rotate.X_AXIS)
    private val rotateY = Rotate(-20.0, Rotate.Y_AXIS)

    // Refresh the unified (collapsed) IDgraph tree.
    private fun refreshTrees() {
        sceneGroup.children.remove(idGraphCollapseGroup)
        // For I tree, internal mapping is identity.
        val (iTreeRender, iVisibleLeaves) = renderCollapsedIDGraph(iTreeRoot, null, globalCollapseLevels, { it.level }, Color.GREEN, 3.0)
        // For D tree, internal mapping: newLevel = n + m - declared level,
        // where n = (max level in I tree) + 1 and m = (max level in D tree) + 1.
        val n = computeMaxLevel(iTreeRoot) + 1
        val m = computeMaxLevel(dTreeRoot) + 1
        val getInternalD: (TreeNode) -> Int = { node -> n + m - node.level }
        val (dTreeRender, dVisibleLeaves) = renderCollapsedIDGraph(dTreeRoot, null, globalCollapseLevels, getInternalD, Color.BLUE, 3.0)
        // Build cross edges (purple) between every visible leaf of I tree and D tree.
        val crossEdges = Group()
        for (iLeaf in iVisibleLeaves) {
            for (dLeaf in dVisibleLeaves) {
                crossEdges.children.add(create3DLine(iLeaf.position, dLeaf.position, 1.0, Color.PURPLE))
            }
        }
        idGraphCollapseGroup = Group(iTreeRender, dTreeRender, crossEdges)
        sceneGroup.children.add(idGraphCollapseGroup)
    }

    // Update the toggle overlay: a single column of toggle arrows for all distinct internal levels.
    private fun updateToggleOverlay() {
        toggleOverlayGroup.children.clear()
        // For I tree, internal levels are 0 .. iMax.
        val iMax = computeMaxLevel(iTreeRoot)
        // For D tree, internal levels (via mapping) range from (n+1) to (n+m) with n = iMax+1.
        val n = iMax + 1
        val m = computeMaxLevel(dTreeRoot) + 1
        val dLevels = (n + 1..n + m).toList()
        val unionLevels = (0..iMax).toList() + dLevels
        val sortedLevels = unionLevels.sorted()
        var i = 0
        for (lvl in sortedLevels) {
            val arrow = if (globalCollapseLevels.contains(lvl)) "▶" else "▼"
            val txt = Text("lvl $lvl: $arrow")
            txt.fill = Color.BLACK
            txt.style = "-fx-font-size: 16px;"
            txt.translateX = -350.0
            txt.translateY = -250.0 + i * 30.0
            txt.setOnMouseClicked {
                if (globalCollapseLevels.contains(lvl)) globalCollapseLevels.remove(lvl)
                else globalCollapseLevels.add(lvl)
                refreshTrees()
                updateToggleOverlay()
            }
            toggleOverlayGroup.children.add(txt)
            i++
        }
    }

    override fun start(primaryStage: Stage) {
        // Define tetrahedron vertices for positioning.
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

        val iTreeDefinition = "A(lvl0, source)->B,C; B(lvl1, int)->D,E; C(lvl1, int)->F; D(lvl2, int)->; E(lvl2,int)->; F(lvl2,int)->"
        val dTreeDefinition = "X(lvl0, source)->Y,Z; Y(lvl1, int)->P; Z(lvl1, int)->Q,R; Q(lvl2,measurement)->; R(lvl2,measurement)->; P(lvl2,measurement)->"

        iTreeRoot = parseTreeDefinition(iTreeDefinition)
        dTreeRoot = parseTreeDefinition(dTreeDefinition)
        assignLeafOrderAndPropagate(iTreeRoot)
        assignLeafOrderAndPropagate(dTreeRoot)
        layoutTree(iTreeRoot, iSource, iLeftEndpoint, iRightEndpoint)
        layoutTree(dTreeRoot, dSource, dLeftEndpoint, dRightEndpoint)

        sceneGroup = Group(edgesGroup, vertexMarkers, vertexLabels)
        refreshTrees()
        sceneGroup.children.add(toggleOverlayGroup)
        updateToggleOverlay()

        val center = Point3D(
            (iSource.x + dSource.x + dLeftEndpoint.x + dRightEndpoint.x) / 4,
            (iSource.y + dSource.y + dLeftEndpoint.y + dRightEndpoint.y) / 4,
            (iSource.z + dSource.z + dLeftEndpoint.z + dRightEndpoint.z) / 4
        )
        sceneGroup.translateX = -center.x
        sceneGroup.translateY = -center.y
        sceneGroup.translateZ = -center.z

        sceneGroup.transforms.addAll(rotateX, rotateY)

        val camera = PerspectiveCamera(true)
        camera.translateZ = -800.0
        camera.nearClip = 0.1
        camera.farClip = 2000.0

        val rootGroup = Group(sceneGroup)
        val scene = Scene(rootGroup, 800.0, 600.0, true, SceneAntialiasing.BALANCED)
        scene.fill = Color.WHITE
        scene.camera = camera

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

        primaryStage.title = "Unified IDgraph Tree with Collapse on All Levels"
        primaryStage.scene = scene
        primaryStage.show()
    }
}

fun main(args: Array<String>) {
    Application.launch(GraphTetrahedronApp::class.java, *args)
}
