// cant collapse
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

// ----------------------------
// Data structures and parser
// ----------------------------

// Tree node with label, level, type, children, a computed position, and a horizontal parameter (u)
data class TreeNode(
    val label: String,
    val level: Int,
    val type: String,
    var children: MutableList<TreeNode> = mutableListOf(),
    var position: Point3D = Point3D(0.0, 0.0, 0.0),
    var u: Double = 0.0  // horizontal parameter for layout (0 = leftmost, 1 = rightmost)
)

// Parse a tree definition string in the format:
// "A(lvl0, source)->B,C; B(lvl1, int)->D; C(lvl1, int)->; D(lvl3, int)"
// For the D‑tree: "X(lvl0, source)->Y,Z; Y(lvl1, int)->P; Z(lvl1, int)->Q,R; Q(lvl2,measurment)->; R(lvl2,measurment)->"
fun parseTreeDefinition(defStr: String): TreeNode {
    val trimmedDef = if (defStr.contains(":")) defStr.substringAfter(":").trim() else defStr.trim()
    if (trimmedDef.isEmpty()) {
        throw IllegalArgumentException("Tree definition string is empty after removing identifier.")
    }
    val segments = trimmedDef.split(";").map { it.trim() }.filter { it.isNotEmpty() }
    if (segments.isEmpty()) {
        throw IllegalArgumentException("No valid segments found in tree definition string.")
    }
    val nodesMap = mutableMapOf<String, TreeNode>()
    val childrenMapping = mutableMapOf<String, List<String>>()
    for (segment in segments) {
        // Example segment: "A(lvl0, source)->B,C"
        val parts = segment.split("->")
        val leftPart = parts[0].trim() // e.g. "A(lvl0, source)"
        val label = leftPart.substringBefore("(").trim()
        val inside = leftPart.substringAfter("(").substringBefore(")").trim() // e.g. "lvl0, source"
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
    // Link children
    for ((parentLabel, childLabels) in childrenMapping) {
        val parentNode = nodesMap[parentLabel]!!
        for (childLabel in childLabels) {
            val childNode = nodesMap.getOrPut(childLabel) { TreeNode(childLabel, 0, "") }
            parentNode.children.add(childNode)
        }
    }
    // Use the first segment's label as the root.
    val rootLabel = segments[0].substringBefore("(").trim()
    return nodesMap[rootLabel]!!
}

// ----------------------------
// Layout Functions
// ----------------------------

// Assign horizontal order (u) to all leaves and propagate upward.
fun assignLeafOrderAndPropagate(root: TreeNode) {
    val leaves = mutableListOf<TreeNode>()
    fun dfs(node: TreeNode) {
        if (node.children.isEmpty()) {
            leaves.add(node)
        } else {
            node.children.forEach { dfs(it) }
        }
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

// Compute the maximum level in the tree.
fun computeMaxLevel(root: TreeNode): Int {
    var maxL = root.level
    for (child in root.children) {
        maxL = maxOf(maxL, computeMaxLevel(child))
    }
    return maxL
}

// Linear interpolation between two Point3D points.
fun interpolate(p1: Point3D, p2: Point3D, t: Double): Point3D {
    return Point3D(
        p1.x + (p2.x - p1.x) * t,
        p1.y + (p2.y - p1.y) * t,
        p1.z + (p2.z - p1.z) * t
    )
}

// Layout the tree on a face defined by a source and two endpoints.
// Each node's position is interpolated based on its level.
fun layoutTree(root: TreeNode, source: Point3D, leftLeaf: Point3D, rightLeaf: Point3D) {
    val maxLevelValue = computeMaxLevel(root).toDouble()
    fun layoutNode(node: TreeNode) {
        val levelFraction = if (maxLevelValue == 0.0) 0.0 else node.level / maxLevelValue
        val baselinePoint = interpolate(leftLeaf, rightLeaf, node.u)
        node.position = interpolate(source, baselinePoint, levelFraction)
        node.children.forEach { layoutNode(it) }
    }
    layoutNode(root)
}

// ----------------------------
// JavaFX Helpers
// ----------------------------

// Create a 3D line (as a Cylinder) between two points.
fun create3DLine(start: Point3D, end: Point3D, radius: Double, color: Color): Cylinder {
    val diff = end.subtract(start)
    val height = diff.magnitude()
    val mid = start.midpoint(end)
    val cylinder = Cylinder(radius, height)
    cylinder.material = PhongMaterial(color)
    // Default alignment is along Y-axis.
    val yAxis = Point3D(0.0, 1.0, 0.0)
    val diffNormalized = diff.normalize()
    val angle = Math.toDegrees(acos(yAxis.dotProduct(diffNormalized)))
    val axisOfRotation = yAxis.crossProduct(diffNormalized)
    if (axisOfRotation.magnitude() > 0) {
        cylinder.transforms.add(Rotate(angle, axisOfRotation))
    }
    cylinder.translateX = mid.x
    cylinder.translateY = mid.y
    cylinder.translateZ = mid.z
    return cylinder
}

// Extension: Compute midpoint between two Point3D objects.
fun Point3D.midpoint(other: Point3D): Point3D {
    return Point3D((this.x + other.x) / 2, (this.y + other.y) / 2, (this.z + other.z) / 2)
}

// Create a sphere marker at a given position.
fun createVertexMarker(position: Point3D, radius: Double, color: Color): Sphere {
    val sphere = Sphere(radius)
    sphere.material = PhongMaterial(color)
    sphere.translateX = position.x
    sphere.translateY = position.y
    sphere.translateZ = position.z
    return sphere
}

// Create a text label at a given position.
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

// Build a JavaFX Group for a tree using each node's actual label.
fun buildTreeGroup(root: TreeNode, nodeRadius: Double, nodeColor: Color): Group {
    val group = Group()
    fun addNode(node: TreeNode) {
        val marker = createVertexMarker(node.position, nodeRadius, nodeColor)
        val label = createTextLabel(node.label, node.position, nodeColor)
        group.children.addAll(marker, label)
        node.children.forEach {
            val edge = create3DLine(node.position, it.position, 1.0, nodeColor)
            group.children.add(edge)
            addNode(it)
        }
    }
    addNode(root)
    return group
}

// ----------------------------
// Helpers to Collect Leaves
// ----------------------------

// Collect leaves that satisfy a given predicate.
fun collectLeavesByType(node: TreeNode, predicate: (TreeNode) -> Boolean, leaves: MutableList<Point3D>) {
    if (node.children.isEmpty()) {
        if (predicate(node)) leaves.add(node.position)
    } else {
        node.children.forEach { collectLeavesByType(it, predicate, leaves) }
    }
}

// Collect all leaves.
fun collectLeaves(node: TreeNode, leaves: MutableList<Point3D>) {
    if (node.children.isEmpty()) {
        leaves.add(node.position)
    } else {
        node.children.forEach { collectLeaves(it, leaves) }
    }
}

// ----------------------------
// Main Application
// ----------------------------

class GraphTetrahedronApp : Application() {

    // For interactive rotation via mouse drag.
    private var anchorX = 0.0
    private var anchorY = 0.0
    private var anchorAngleX = 0.0
    private var anchorAngleY = 0.0

    // Two rotation transforms applied to the 3D scene.
    private val rotateX = Rotate(20.0, Rotate.X_AXIS)
    private val rotateY = Rotate(-20.0, Rotate.Y_AXIS)

    override fun start(primaryStage: Stage) {
        // --- Define Tetrahedron Vertices ---
        // These vertices are used for the tetrahedron and to position the trees.
        val iSource = Point3D(0.0, 100.0, 0.0)            // Source for I tree
        val dSource = Point3D(-100.0, -100.0, 100.0)         // Source for D tree
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

        // --- Build Tetrahedron Edges and Vertex Markers ---
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

        // --- Define Tree Specification Strings ---
        // I-tree: nodes A, B, C, D
        val iTreeDefinition = "A(lvl0, source)->B,C; B(lvl1, int)->D,E; C(lvl1, int)->F; D(lvl2, int)->;E(lvl2,int)->;F(lvl2,int)->"
        // D-tree: nodes X, Y, Z, P, Q, R
        val dTreeDefinition = "X(lvl0, source)->Y,Z; Y(lvl1, int)->P; Z(lvl1, int)->Q,R; Q(lvl2,measurment)->; R(lvl2,measurment)->; P(lvl2, measurment)->"

        // --- Parse Tree Definitions ---
        val iTreeRoot = parseTreeDefinition(iTreeDefinition)
        val dTreeRoot = parseTreeDefinition(dTreeDefinition)

        // --- Assign Horizontal Order and Layout ---
        assignLeafOrderAndPropagate(iTreeRoot)
        assignLeafOrderAndPropagate(dTreeRoot)
        layoutTree(iTreeRoot, iSource, iLeftEndpoint, iRightEndpoint)
        layoutTree(dTreeRoot, dSource, dLeftEndpoint, dRightEndpoint)

        // --- Build JavaFX Groups for the Trees (using actual node labels) ---
        val iTreeGroup = buildTreeGroup(iTreeRoot, 3.0, Color.GREEN)
        val dTreeGroup = buildTreeGroup(dTreeRoot, 3.0, Color.BLUE)

        // --- Connect I-tree Leaves to D-tree Measurement Nodes ---
        val iLeaves = mutableListOf<Point3D>()
        collectLeaves(iTreeRoot, iLeaves)
        val dMeasurementLeaves = mutableListOf<Point3D>()
        collectLeavesByType(dTreeRoot, { it.type.equals("measurment", ignoreCase = true) }, dMeasurementLeaves)
        val connectionLines = Group()
        for (iLeaf in iLeaves) {
            for (dLeaf in dMeasurementLeaves) {
                val connection = create3DLine(iLeaf, dLeaf, 0.5, Color.GREEN)
                connection.viewOrder = -1.0
                connectionLines.children.add(connection)
            }
        }

        // --- Center the Entire Tetrahedron ---
        val center = Point3D(
            (iSource.x + dSource.x + dLeftEndpoint.x + dRightEndpoint.x) / 4,
            (iSource.y + dSource.y + dLeftEndpoint.y + dRightEndpoint.y) / 4,
            (iSource.z + dSource.z + dLeftEndpoint.z + dRightEndpoint.z) / 4
        )
        val sceneGroup = Group(edgesGroup, vertexMarkers, vertexLabels, dTreeGroup, iTreeGroup, connectionLines)
        sceneGroup.translateX = -center.x
        sceneGroup.translateY = -center.y
        sceneGroup.translateZ = -center.z
        sceneGroup.transforms.addAll(rotateX, rotateY)

        // --- Set Up Camera and Scene ---
        val camera = PerspectiveCamera(true)
        camera.translateZ = -800.0
        camera.nearClip = 0.1
        camera.farClip = 2000.0

        val rootGroup = Group(sceneGroup)
        val scene = Scene(rootGroup, 800.0, 600.0, true, SceneAntialiasing.BALANCED)
        scene.fill = Color.WHITE
        scene.camera = camera

        // --- Interactivity: Mouse Drag for Rotation, Scroll for Zoom ---
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

        // --- Animation Timer to Keep Labels Facing the Camera ---
        object : AnimationTimer() {
            override fun handle(now: Long) {
                for (node in sceneGroup.children) {
                    if (node is Text) {
                        node.transforms.setAll(
                            Rotate(-rotateY.angle, Rotate.Y_AXIS),
                            Rotate(-rotateX.angle, Rotate.X_AXIS)
                        )
                    }
                }
            }
        }.start()

        primaryStage.title = "Generalized Tetrahedron with Custom Tree Definitions"
        primaryStage.scene = scene
        primaryStage.show()
    }
}

fun main(args: Array<String>) {
    Application.launch(GraphTetrahedronApp::class.java, *args)
}


// 2 can collapse but not measurment + errors when collapsing adjescent lvls
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

// ----------------------------
// Data Structures and Parser
// ----------------------------

// Tree node with label, level, type, children, a computed position, and a horizontal parameter (u)
data class TreeNode(
    val label: String,
    val level: Int,
    val type: String,
    var children: MutableList<TreeNode> = mutableListOf(),
    var position: Point3D = Point3D(0.0, 0.0, 0.0),
    var u: Double = 0.0  // horizontal parameter for layout (0 = leftmost, 1 = rightmost)
)

// Parse a tree definition string in the format:
// "A(lvl0, source)->B,C; B(lvl1, int)->D; C(lvl1, int)->; D(lvl3, int)"
fun parseTreeDefinition(defStr: String): TreeNode {
    val trimmedDef = if (defStr.contains(":")) defStr.substringAfter(":").trim() else defStr.trim()
    if (trimmedDef.isEmpty()) {
        throw IllegalArgumentException("Tree definition string is empty after removing identifier.")
    }
    val segments = trimmedDef.split(";").map { it.trim() }.filter { it.isNotEmpty() }
    if (segments.isEmpty()) {
        throw IllegalArgumentException("No valid segments found in tree definition string.")
    }
    val nodesMap = mutableMapOf<String, TreeNode>()
    val childrenMapping = mutableMapOf<String, List<String>>()
    for (segment in segments) {
        // Example: "A(lvl0, source)->B,C"
        val parts = segment.split("->")
        val leftPart = parts[0].trim() // e.g. "A(lvl0, source)"
        val label = leftPart.substringBefore("(").trim()
        val inside = leftPart.substringAfter("(").substringBefore(")").trim() // e.g. "lvl0, source"
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
    // Link children
    for ((parentLabel, childLabels) in childrenMapping) {
        val parentNode = nodesMap[parentLabel]!!
        for (childLabel in childLabels) {
            val childNode = nodesMap.getOrPut(childLabel) { TreeNode(childLabel, 0, "") }
            parentNode.children.add(childNode)
        }
    }
    // Use the first segment's label as the root.
    val rootLabel = segments[0].substringBefore("(").trim()
    return nodesMap[rootLabel]!!
}

// ----------------------------
// Layout Functions
// ----------------------------

// Assign horizontal order (u) to all leaves and propagate upward.
fun assignLeafOrderAndPropagate(root: TreeNode) {
    val leaves = mutableListOf<TreeNode>()
    fun dfs(node: TreeNode) {
        if (node.children.isEmpty()) {
            leaves.add(node)
        } else {
            node.children.forEach { dfs(it) }
        }
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

// Compute the maximum level in the tree.
fun computeMaxLevel(root: TreeNode): Int {
    var maxL = root.level
    for (child in root.children) {
        maxL = maxOf(maxL, computeMaxLevel(child))
    }
    return maxL
}

// Linear interpolation between two Point3D points.
fun interpolate(p1: Point3D, p2: Point3D, t: Double): Point3D {
    return Point3D(
        p1.x + (p2.x - p1.x) * t,
        p1.y + (p2.y - p1.y) * t,
        p1.z + (p2.z - p1.z) * t
    )
}

// Layout the tree on a face defined by a source and two endpoints.
fun layoutTree(root: TreeNode, source: Point3D, leftLeaf: Point3D, rightLeaf: Point3D) {
    val maxLevelValue = computeMaxLevel(root).toDouble()
    fun layoutNode(node: TreeNode) {
        val levelFraction = if (maxLevelValue == 0.0) 0.0 else node.level / maxLevelValue
        val baselinePoint = interpolate(leftLeaf, rightLeaf, node.u)
        node.position = interpolate(source, baselinePoint, levelFraction)
        node.children.forEach { layoutNode(it) }
    }
    layoutNode(root)
}

// ----------------------------
// JavaFX Helpers
// ----------------------------

// Create a 3D line (as a Cylinder) between two points.
fun create3DLine(start: Point3D, end: Point3D, radius: Double, color: Color): Cylinder {
    val diff = end.subtract(start)
    val height = diff.magnitude()
    val mid = start.midpoint(end)
    val cylinder = Cylinder(radius, height)
    cylinder.material = PhongMaterial(color)
    val yAxis = Point3D(0.0, 1.0, 0.0)
    val diffNormalized = diff.normalize()
    val angle = Math.toDegrees(acos(yAxis.dotProduct(diffNormalized)))
    val axisOfRotation = yAxis.crossProduct(diffNormalized)
    if (axisOfRotation.magnitude() > 0) {
        cylinder.transforms.add(Rotate(angle, axisOfRotation))
    }
    cylinder.translateX = mid.x
    cylinder.translateY = mid.y
    cylinder.translateZ = mid.z
    return cylinder
}

// Extension: Compute midpoint between two Point3D objects.
fun Point3D.midpoint(other: Point3D): Point3D {
    return Point3D((this.x + other.x) / 2, (this.y + other.y) / 2, (this.z + other.z) / 2)
}

// Create a sphere marker at a given position.
fun createVertexMarker(position: Point3D, radius: Double, color: Color): Sphere {
    val sphere = Sphere(radius)
    sphere.material = PhongMaterial(color)
    sphere.translateX = position.x
    sphere.translateY = position.y
    sphere.translateZ = position.z
    return sphere
}

// Create a text label at a given position.
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

// ----------------------------
// Rendering with Level Collapse
// ----------------------------
// This recursive function renders a tree node while “collapsing” any node whose level is in the
// provided collapsedLevels set. When a node is collapsed its marker/label is not rendered; instead,
// shortcut edges are drawn from its effective parent (the nearest rendered node from the level above)
// directly to each of its children.
fun renderNodeWithCollapses(
    node: TreeNode,
    effectiveParent: TreeNode?,
    collapsedLevels: Set<Int>,
    nodeRadius: Double,
    nodeColor: Color
): Group {
    val group = Group()
    if (node.level in collapsedLevels) {
        // Node is collapsed – do not render its marker/label.
        // For each child, if an effective parent exists, draw an edge from it to the child.
        for (child in node.children) {
            if (effectiveParent != null) {
                val edge = create3DLine(effectiveParent.position, child.position, 1.0, nodeColor)
                group.children.add(edge)
            }
            group.children.add(renderNodeWithCollapses(child, effectiveParent, collapsedLevels, nodeRadius, nodeColor))
        }
    } else {
        // Render the node normally.
        val marker = createVertexMarker(node.position, nodeRadius, nodeColor)
        val label = createTextLabel(node.label, node.position, nodeColor)
        group.children.addAll(marker, label)
        for (child in node.children) {
            if (child.level in collapsedLevels) {
                // Instead of rendering the collapsed child, draw shortcut edges from the current node
                // to each grandchild.
                for (grandchild in child.children) {
                    val edge = create3DLine(node.position, grandchild.position, 1.0, nodeColor)
                    group.children.add(edge)
                    group.children.add(renderNodeWithCollapses(grandchild, node, collapsedLevels, nodeRadius, nodeColor))
                }
            } else {
                val edge = create3DLine(node.position, child.position, 1.0, nodeColor)
                group.children.add(edge)
                group.children.add(renderNodeWithCollapses(child, node, collapsedLevels, nodeRadius, nodeColor))
            }
        }
    }
    return group
}

// External function to render an entire tree with the given set of collapsed levels.
fun renderTreeWithCollapse(tree: TreeNode, collapsedLevels: Set<Int>, nodeRadius: Double, nodeColor: Color): Group {
    return renderNodeWithCollapses(tree, null, collapsedLevels, nodeRadius, nodeColor)
}

// ----------------------------
// Main Application
// ----------------------------
class GraphTetrahedronApp : Application() {

    // The tree roots (unchanged) and scene groups.
    lateinit var iTreeRoot: TreeNode
    lateinit var dTreeRoot: TreeNode
    lateinit var sceneGroup: Group
    lateinit var overlayGroup: Group
    var iTreeGroup: Group = Group()
    var dTreeGroup: Group = Group()

    // Sets to track which levels are collapsed (per tree)
    val collapsedLevelsI: MutableSet<Int> = mutableSetOf()
    val collapsedLevelsD: MutableSet<Int> = mutableSetOf()

    // Groups for the toggle controls (displayed as overlay text)
    var iToggleGroup: Group = Group()
    var dToggleGroup: Group = Group()

    // For interactive rotation via mouse drag.
    private var anchorX = 0.0
    private var anchorY = 0.0
    private var anchorAngleX = 0.0
    private var anchorAngleY = 0.0

    // Two rotation transforms applied to the 3D scene.
    private val rotateX = Rotate(20.0, Rotate.X_AXIS)
    private val rotateY = Rotate(-20.0, Rotate.Y_AXIS)

    // Refresh the tree renderings.
    fun refreshTrees() {
        sceneGroup.children.remove(iTreeGroup)
        sceneGroup.children.remove(dTreeGroup)
        iTreeGroup = renderTreeWithCollapse(iTreeRoot, collapsedLevelsI, 3.0, Color.GREEN)
        dTreeGroup = renderTreeWithCollapse(dTreeRoot, collapsedLevelsD, 3.0, Color.BLUE)
        sceneGroup.children.addAll(iTreeGroup, dTreeGroup)
    }

    // Update the toggle overlay controls.
    fun updateToggleOverlay() {
        iToggleGroup.children.clear()
        dToggleGroup.children.clear()

        val iMaxLevel = computeMaxLevel(iTreeRoot)
        val dMaxLevel = computeMaxLevel(dTreeRoot)

        // For I tree: display toggles on the left.
        for (level in 0..iMaxLevel) {
            val isCollapsed = level in collapsedLevelsI
            val arrow = if (isCollapsed) "▶" else "▼"
            val toggleText = Text("I lvl$level: $arrow")
            toggleText.fill = Color.BLACK
            toggleText.style = "-fx-font-size: 16px;"
            toggleText.translateX = -350.0
            toggleText.translateY = -250.0 + level * 30
            toggleText.setOnMouseClicked {
                if (collapsedLevelsI.contains(level)) {
                    collapsedLevelsI.remove(level)
                } else {
                    collapsedLevelsI.add(level)
                }
                refreshTrees()
                updateToggleOverlay()
            }
            iToggleGroup.children.add(toggleText)
        }

        // For D tree: display toggles on the right.
        for (level in 0..dMaxLevel) {
            val isCollapsed = level in collapsedLevelsD
            val arrow = if (isCollapsed) "▶" else "▼"
            val toggleText = Text("D lvl$level: $arrow")
            toggleText.fill = Color.BLACK
            toggleText.style = "-fx-font-size: 16px;"
            toggleText.translateX = 250.0
            toggleText.translateY = -250.0 + level * 30
            toggleText.setOnMouseClicked {
                if (collapsedLevelsD.contains(level)) {
                    collapsedLevelsD.remove(level)
                } else {
                    collapsedLevelsD.add(level)
                }
                refreshTrees()
                updateToggleOverlay()
            }
            dToggleGroup.children.add(toggleText)
        }
    }

    override fun start(primaryStage: Stage) {
        // --- Define Tetrahedron Vertices ---
        val iSource = Point3D(0.0, 100.0, 0.0)
        val dSource = Point3D(-100.0, -100.0, 100.0)
        val dLeftEndpoint = Point3D(100.0, -100.0, 100.0)
        val dRightEndpoint = Point3D(0.0, -100.0, -100.0)
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

        // --- Build Tetrahedron Edges and Vertex Markers ---
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

        // --- Define Tree Specification Strings (unchanged) ---
        val iTreeDefinition = "A(lvl0, source)->B,C; B(lvl1, int)->D,E; C(lvl1, int)->F; D(lvl2, int)->;E(lvl2,int)->;F(lvl2,int)->"
        val dTreeDefinition = "X(lvl0, source)->Y,Z; Y(lvl1, int)->P; Z(lvl1, int)->Q,R; Q(lvl2,int)->M,N; R(lvl3,measurment)->;P(lvl3, measurment)->;M(lvl3, measurment)->;N(lvl3,measurment)->"

        // --- Parse Tree Definitions ---
        iTreeRoot = parseTreeDefinition(iTreeDefinition)
        dTreeRoot = parseTreeDefinition(dTreeDefinition)

        // --- Assign Horizontal Order and Layout ---
        assignLeafOrderAndPropagate(iTreeRoot)
        assignLeafOrderAndPropagate(dTreeRoot)
        layoutTree(iTreeRoot, iSource, iLeftEndpoint, iRightEndpoint)
        layoutTree(dTreeRoot, dSource, dLeftEndpoint, dRightEndpoint)

        // --- Build Initial Tree Groups using renderTreeWithCollapse ---
        iTreeGroup = renderTreeWithCollapse(iTreeRoot, collapsedLevelsI, 3.0, Color.GREEN)
        dTreeGroup = renderTreeWithCollapse(dTreeRoot, collapsedLevelsD, 3.0, Color.BLUE)

        // --- Connect I-tree Leaves to D-tree Measurement Nodes (unchanged) ---
        val iLeaves = mutableListOf<Point3D>()
        fun collectLeaves(node: TreeNode, leaves: MutableList<Point3D>) {
            if (node.children.isEmpty()) leaves.add(node.position)
            else node.children.forEach { collectLeaves(it, leaves) }
        }
        collectLeaves(iTreeRoot, iLeaves)
        val dMeasurementLeaves = mutableListOf<Point3D>()
        fun collectLeavesByType(node: TreeNode, predicate: (TreeNode) -> Boolean, leaves: MutableList<Point3D>) {
            if (node.children.isEmpty()) {
                if (predicate(node)) leaves.add(node.position)
            } else {
                node.children.forEach { collectLeavesByType(it, predicate, leaves) }
            }
        }
        collectLeavesByType(dTreeRoot, { it.type.equals("measurment", ignoreCase = true) }, dMeasurementLeaves)
        val connectionLines = Group()
        for (iLeaf in iLeaves) {
            for (dLeaf in dMeasurementLeaves) {
                val connection = create3DLine(iLeaf, dLeaf, 0.5, Color.GREEN)
                connection.viewOrder = -1.0
                connectionLines.children.add(connection)
            }
        }

        // --- Build Scene Group ---
        sceneGroup = Group(edgesGroup, vertexMarkers, vertexLabels, iTreeGroup, dTreeGroup, connectionLines)
        // Center the tetrahedron
        val center = Point3D(
            (iSource.x + dSource.x + dLeftEndpoint.x + dRightEndpoint.x) / 4,
            (iSource.y + dSource.y + dLeftEndpoint.y + dRightEndpoint.y) / 4,
            (iSource.z + dSource.z + dLeftEndpoint.z + dRightEndpoint.z) / 4
        )
        sceneGroup.translateX = -center.x
        sceneGroup.translateY = -center.y
        sceneGroup.translateZ = -center.z
        sceneGroup.transforms.addAll(rotateX, rotateY)

        // --- Create Overlay for Toggle Controls ---
        overlayGroup = Group()
        iToggleGroup = Group()
        dToggleGroup = Group()
        updateToggleOverlay() // Build toggle labels for both trees
        overlayGroup.children.addAll(iToggleGroup, dToggleGroup)

        // --- Build Root Group ---
        val rootGroup = Group(sceneGroup, overlayGroup)

        // --- Set Up Camera and Scene ---
        val camera = PerspectiveCamera(true)
        camera.translateZ = -800.0
        camera.nearClip = 0.1
        camera.farClip = 2000.0

        val scene = Scene(rootGroup, 800.0, 600.0, true, SceneAntialiasing.BALANCED)
        scene.fill = Color.WHITE
        scene.camera = camera

        // --- Interactivity: Mouse Drag for Rotation, Scroll for Zoom ---
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

        // --- Animation Timer to Keep 3D Labels Facing the Camera ---
        object : AnimationTimer() {
            override fun handle(now: Long) {
                for (node in sceneGroup.children) {
                    if (node is Text) {
                        node.transforms.setAll(
                            Rotate(-rotateY.angle, Rotate.Y_AXIS),
                            Rotate(-rotateX.angle, Rotate.X_AXIS)
                        )
                    }
                }
            }
        }.start()

        primaryStage.title = "Generalized Tetrahedron with Level Collapse Toggles"
        primaryStage.scene = scene
        primaryStage.show()
    }
}

fun main(args: Array<String>) {
    Application.launch(GraphTetrahedronApp::class.java, *args)
}

// 3. can collapse adj lvls properly, still problem with measurment lvl
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

// ----------------------------
// Data Structures and Parser
// ----------------------------
data class TreeNode(
    val label: String,
    val level: Int,
    val type: String,
    var children: MutableList<TreeNode> = mutableListOf(),
    var position: Point3D = Point3D(0.0, 0.0, 0.0),
    var u: Double = 0.0 // horizontal parameter (0 = left, 1 = right)
)

fun parseTreeDefinition(defStr: String): TreeNode {
    // Remove any prefix before the colon if present.
    val trimmedDef = if (defStr.contains(":")) defStr.substringAfter(":").trim() else defStr.trim()
    if (trimmedDef.isEmpty()) throw IllegalArgumentException("Empty definition!")
    val segments = trimmedDef.split(";").map { it.trim() }.filter { it.isNotEmpty() }
    val nodesMap = mutableMapOf<String, TreeNode>()
    val childrenMapping = mutableMapOf<String, List<String>>()
    for (segment in segments) {
        // Example segment: "A(lvl0, source)->B,C"
        val parts = segment.split("->")
        val leftPart = parts[0].trim()
        val label = leftPart.substringBefore("(").trim()
        val inside = leftPart.substringAfter("(").substringBefore(")").trim()
        val tokens = inside.split(",").map { it.trim() }
        val level = tokens[0].substringAfter("lvl").toIntOrNull() ?: 0
        val type = if (tokens.size > 1) tokens[1] else ""
        nodesMap[label] = TreeNode(label, level, type)
        if (parts.size > 1) {
            val childLabels = parts[1].split(",").map { it.trim() }.filter { it.isNotEmpty() }
            childrenMapping[label] = childLabels
        }
    }
    // Link children
    for ((parent, childLabels) in childrenMapping) {
        val parentNode = nodesMap[parent]!!
        for (child in childLabels) {
            val childNode = nodesMap.getOrPut(child) { TreeNode(child, 0, "") }
            parentNode.children.add(childNode)
        }
    }
    // Return the root (first segment)
    val rootLabel = segments[0].substringBefore("(").trim()
    return nodesMap[rootLabel]!!
}

// ----------------------------
// Layout Functions
// ----------------------------
fun assignLeafOrderAndPropagate(root: TreeNode) {
    val leaves = mutableListOf<TreeNode>()
    fun dfs(node: TreeNode) {
        if (node.children.isEmpty()) leaves.add(node)
        else node.children.forEach { dfs(it) }
    }
    dfs(root)
    val n = leaves.size
    leaves.forEachIndexed { index, leaf ->
        leaf.u = if (n > 1) index.toDouble() / (n - 1) else 0.5
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

fun interpolate(p1: Point3D, p2: Point3D, t: Double): Point3D =
    Point3D(
        p1.x + (p2.x - p1.x) * t,
        p1.y + (p2.y - p1.y) * t,
        p1.z + (p2.z - p1.z) * t
    )

fun layoutTree(root: TreeNode, source: Point3D, leftLeaf: Point3D, rightLeaf: Point3D) {
    val maxLevelVal = computeMaxLevel(root).toDouble()
    fun layout(node: TreeNode) {
        val t = if (maxLevelVal == 0.0) 0.0 else node.level / maxLevelVal
        val baseline = interpolate(leftLeaf, rightLeaf, node.u)
        node.position = interpolate(source, baseline, t)
        node.children.forEach { layout(it) }
    }
    layout(root)
}

// ----------------------------
// JavaFX 3D Helpers
// ----------------------------
fun create3DLine(start: Point3D, end: Point3D, radius: Double, color: Color): Cylinder {
    val diff = end.subtract(start)
    val height = diff.magnitude()
    val mid = Point3D((start.x + end.x) / 2, (start.y + end.y) / 2, (start.z + end.z) / 2)
    val cyl = Cylinder(radius, height)
    cyl.material = PhongMaterial(color)
    val yAxis = Point3D(0.0, 1.0, 0.0)
    val diffNorm = diff.normalize()
    val angle = Math.toDegrees(acos(yAxis.dotProduct(diffNorm)))
    val axis = yAxis.crossProduct(diffNorm)
    if (axis.magnitude() > 0)
        cyl.transforms.add(Rotate(angle, axis))
    cyl.translateX = mid.x
    cyl.translateY = mid.y
    cyl.translateZ = mid.z
    return cyl
}

fun createVertexMarker(pos: Point3D, radius: Double, color: Color): Sphere {
    val s = Sphere(radius)
    s.material = PhongMaterial(color)
    s.translateX = pos.x
    s.translateY = pos.y
    s.translateZ = pos.z
    return s
}

fun createTextLabel(textStr: String, pos: Point3D, color: Color): Text {
    val t = Text(textStr)
    t.fill = color
    t.translateX = pos.x
    t.translateY = pos.y
    t.translateZ = pos.z
    t.scaleX = 2.0
    t.scaleY = 2.0
    return t
}

// ----------------------------
// Revised Collapse Rendering
// ----------------------------
/*
  The function below recursively renders the tree. Each node is checked for “visibility.”
  A node is considered collapsed if its level is in the collapse set OR (if its type is
  "measurment") when measurementCollapse is true. In that case, the node's marker and label
  are skipped and its children are rendered as if attached to the last visible (effective) parent.
*/
fun renderTree(
    node: TreeNode,
    effectiveParent: TreeNode?,
    collapseLevels: Set<Int>,
    measurementCollapse: Boolean,
    nodeRadius: Double,
    nodeColor: Color
): Group {
    val group = Group()
    val isCollapsed = (node.level in collapseLevels) ||
            (node.type.equals("measurment", ignoreCase = true) && measurementCollapse)
    // If the current node is visible, draw it.
    if (!isCollapsed) {
        group.children.add(createVertexMarker(node.position, nodeRadius, nodeColor))
        group.children.add(createTextLabel(node.label, node.position, nodeColor))
    }
    // The effective parent is the last visible node.
    val newEffective = if (!isCollapsed) node else effectiveParent

    // Process children: if a child is visible, draw an edge from the effective parent.
    for (child in node.children) {
        val childCollapsed = (child.level in collapseLevels) ||
                (child.type.equals("measurment", ignoreCase = true) && measurementCollapse)
        if (newEffective != null && !childCollapsed) {
            group.children.add(create3DLine(newEffective.position, child.position, 1.0, nodeColor))
        }
        group.children.add(renderTree(child, newEffective, collapseLevels, measurementCollapse, nodeRadius, nodeColor))
    }
    return group
}

fun renderTreeWithCollapse(
    root: TreeNode,
    collapseLevels: Set<Int>,
    measurementCollapse: Boolean,
    nodeRadius: Double,
    nodeColor: Color
): Group =
    renderTree(root, null, collapseLevels, measurementCollapse, nodeRadius, nodeColor)

// ----------------------------
// Main Application
// ----------------------------
class GraphTetrahedronApp : Application() {

    // Trees and their root nodes.
    lateinit var iTreeRoot: TreeNode
    lateinit var dTreeRoot: TreeNode

    // Scene groups.
    lateinit var sceneGroup: Group
    lateinit var overlayGroup: Group
    var iTreeGroup: Group = Group()
    var dTreeGroup: Group = Group()

    // Sets to track which levels are collapsed.
    val collapsedLevelsI: MutableSet<Int> = mutableSetOf()
    val collapsedLevelsD: MutableSet<Int> = mutableSetOf()

    // Toggle for D tree "measurment" nodes.
    var dMeasurementCollapsed: Boolean = false

    // Overlay groups for toggles.
    var iToggleGroup: Group = Group()
    var dToggleGroup: Group = Group()

    // For mouse rotation.
    private var anchorX = 0.0
    private var anchorY = 0.0
    private var anchorAngleX = 0.0
    private var anchorAngleY = 0.0
    private val rotateX = Rotate(20.0, Rotate.X_AXIS)
    private val rotateY = Rotate(-20.0, Rotate.Y_AXIS)

    fun refreshTrees() {
        sceneGroup.children.remove(iTreeGroup)
        sceneGroup.children.remove(dTreeGroup)
        iTreeGroup = renderTreeWithCollapse(iTreeRoot, collapsedLevelsI, false, 3.0, Color.GREEN)
        dTreeGroup = renderTreeWithCollapse(dTreeRoot, collapsedLevelsD, dMeasurementCollapsed, 3.0, Color.BLUE)
        sceneGroup.children.addAll(iTreeGroup, dTreeGroup)
    }

    fun updateToggleOverlay() {
        iToggleGroup.children.clear()
        dToggleGroup.children.clear()

        val iMax = computeMaxLevel(iTreeRoot)
        val dMax = computeMaxLevel(dTreeRoot)

        // I-tree toggles (left side)
        for (lvl in 0..iMax) {
            val arrow = if (lvl in collapsedLevelsI) "▶" else "▼"
            val txt = Text("I lvl$lvl: $arrow")
            txt.fill = Color.BLACK
            txt.style = "-fx-font-size: 16px;"
            txt.translateX = -350.0
            txt.translateY = -250.0 + lvl * 30
            txt.setOnMouseClicked {
                if (lvl in collapsedLevelsI) collapsedLevelsI.remove(lvl)
                else collapsedLevelsI.add(lvl)
                refreshTrees()
                updateToggleOverlay()
            }
            iToggleGroup.children.add(txt)
        }

        // D-tree level toggles (right side)
        for (lvl in 0..dMax) {
            val arrow = if (lvl in collapsedLevelsD) "▶" else "▼"
            val txt = Text("D lvl$lvl: $arrow")
            txt.fill = Color.BLACK
            txt.style = "-fx-font-size: 16px;"
            txt.translateX = 250.0
            txt.translateY = -250.0 + lvl * 30
            txt.setOnMouseClicked {
                if (lvl in collapsedLevelsD) collapsedLevelsD.remove(lvl)
                else collapsedLevelsD.add(lvl)
                refreshTrees()
                updateToggleOverlay()
            }
            dToggleGroup.children.add(txt)
        }
        // Extra toggle for measurement nodes in D-tree.
        val measArrow = if (dMeasurementCollapsed) "▶" else "▼"
        val measTxt = Text("D Measurement: $measArrow")
        measTxt.fill = Color.BLACK
        measTxt.style = "-fx-font-size: 16px;"
        measTxt.translateX = 250.0
        measTxt.translateY = -250.0 + (dMax + 1) * 30
        measTxt.setOnMouseClicked {
            dMeasurementCollapsed = !dMeasurementCollapsed
            refreshTrees()
            updateToggleOverlay()
        }
        dToggleGroup.children.add(measTxt)
    }

    override fun start(primaryStage: Stage) {
        // Define tetrahedron vertices.
        val iSource = Point3D(0.0, 100.0, 0.0)
        val dSource = Point3D(-100.0, -100.0, 100.0)
        val dLeft = Point3D(100.0, -100.0, 100.0)
        val dRight = Point3D(0.0, -100.0, -100.0)
        val iLeft = Point3D(
            iSource.x + 0.5 * (dLeft.x - iSource.x),
            iSource.y + 0.5 * (dLeft.y - iSource.y),
            iSource.z + 0.5 * (dLeft.z - iSource.z)
        )
        val iRight = Point3D(
            iSource.x + 0.5 * (dRight.x - iSource.x),
            iSource.y + 0.5 * (dRight.y - iSource.y),
            iSource.z + 0.5 * (dRight.z - iSource.z)
        )

        // Build tetrahedron edges and markers.
        val tetraColor = Color.LIGHTGRAY.deriveColor(0.0, 1.0, 1.0, 0.3)
        val edge1 = create3DLine(iSource, dSource, 1.0, tetraColor)
        val edge2 = create3DLine(iSource, dLeft, 1.0, tetraColor)
        val edge3 = create3DLine(iSource, dRight, 1.0, tetraColor)
        val edge4 = create3DLine(dSource, dLeft, 1.0, tetraColor)
        val edge5 = create3DLine(dSource, dRight, 1.0, tetraColor)
        val edge6 = create3DLine(dLeft, dRight, 1.0, tetraColor)
        val edgesGroup = Group(edge1, edge2, edge3, edge4, edge5, edge6)
        edgesGroup.viewOrder = 0.0

        val vertexMarkers = Group(
            createVertexMarker(iSource, 5.0, Color.GREEN),
            createVertexMarker(dSource, 5.0, Color.ORANGE),
            createVertexMarker(dLeft, 5.0, Color.ORANGE),
            createVertexMarker(dRight, 5.0, Color.ORANGE)
        )
        val vertexLabels = Group(
            createTextLabel("I", iSource, Color.GREEN),
            createTextLabel("d", dSource, Color.ORANGE),
            createTextLabel("m1", dLeft, Color.ORANGE),
            createTextLabel("m2", dRight, Color.ORANGE)
        )

        // Define tree definitions.
        val iTreeDef = "A(lvl0, source)->B,C; B(lvl1, int)->D,E; C(lvl1, int)->F; D(lvl2, int)->; E(lvl2, int)->; F(lvl2, int)->"
        val dTreeDef = "X(lvl0, source)->Y,Z; Y(lvl1, int)->P; Z(lvl1, int)->Q,R; Q(lvl2, int)->M,N; R(lvl3, measurment)->; P(lvl3, measurment)->; M(lvl3, measurment)->; N(lvl3, measurment)->"

        // Parse trees.
        iTreeRoot = parseTreeDefinition(iTreeDef)
        dTreeRoot = parseTreeDefinition(dTreeDef)

        // Compute horizontal orders and layout.
        assignLeafOrderAndPropagate(iTreeRoot)
        assignLeafOrderAndPropagate(dTreeRoot)
        layoutTree(iTreeRoot, iSource, iLeft, iRight)
        layoutTree(dTreeRoot, dSource, dLeft, dRight)

        // Build tree groups.
        iTreeGroup = renderTreeWithCollapse(iTreeRoot, collapsedLevelsI, false, 3.0, Color.GREEN)
        dTreeGroup = renderTreeWithCollapse(dTreeRoot, collapsedLevelsD, dMeasurementCollapsed, 3.0, Color.BLUE)

        // Optionally, add connections between I-tree leaves and D-tree measurement nodes.
        val iLeaves = mutableListOf<Point3D>()
        fun collectLeaves(node: TreeNode, list: MutableList<Point3D>) {
            if (node.children.isEmpty()) list.add(node.position)
            else node.children.forEach { collectLeaves(it, list) }
        }
        collectLeaves(iTreeRoot, iLeaves)
        val dMeasLeaves = mutableListOf<Point3D>()
        fun collectMeasurementLeaves(node: TreeNode, list: MutableList<Point3D>) {
            if (node.children.isEmpty() && node.type.equals("measurment", ignoreCase = true))
                list.add(node.position)
            else node.children.forEach { collectMeasurementLeaves(it, list) }
        }
        collectMeasurementLeaves(dTreeRoot, dMeasLeaves)
        val connections = Group()
        for (iLeaf in iLeaves)
            for (dLeaf in dMeasLeaves) {
                val line = create3DLine(iLeaf, dLeaf, 0.5, Color.GREEN)
                line.viewOrder = -1.0
                connections.children.add(line)
            }

        // Build the main scene group.
        sceneGroup = Group(edgesGroup, vertexMarkers, vertexLabels, iTreeGroup, dTreeGroup, connections)
        // Center the tetrahedron.
        val center = Point3D(
            (iSource.x + dSource.x + dLeft.x + dRight.x) / 4,
            (iSource.y + dSource.y + dLeft.y + dRight.y) / 4,
            (iSource.z + dSource.z + dLeft.z + dRight.z) / 4
        )
        sceneGroup.translateX = -center.x
        sceneGroup.translateY = -center.y
        sceneGroup.translateZ = -center.z
        sceneGroup.transforms.addAll(rotateX, rotateY)

        // Create overlay for toggles.
        overlayGroup = Group()
        iToggleGroup = Group()
        dToggleGroup = Group()
        updateToggleOverlay()
        overlayGroup.children.addAll(iToggleGroup, dToggleGroup)

        // Root group.
        val rootGroup = Group(sceneGroup, overlayGroup)

        // Setup camera and scene.
        val camera = PerspectiveCamera(true)
        camera.translateZ = -800.0
        camera.nearClip = 0.1
        camera.farClip = 2000.0
        val scene = Scene(rootGroup, 800.0, 600.0, true, SceneAntialiasing.BALANCED)
        scene.fill = Color.WHITE
        scene.camera = camera

        // Mouse interactivity for rotation and zoom.
        scene.addEventHandler(MouseEvent.MOUSE_PRESSED) { e ->
            anchorX = e.sceneX
            anchorY = e.sceneY
            anchorAngleX = rotateX.angle
            anchorAngleY = rotateY.angle
        }
        scene.addEventHandler(MouseEvent.MOUSE_DRAGGED) { e ->
            rotateX.angle = anchorAngleX - (e.sceneY - anchorY)
            rotateY.angle = anchorAngleY + (e.sceneX - anchorX)
        }
        scene.addEventHandler(ScrollEvent.SCROLL) { e ->
            camera.translateZ += e.deltaY
        }

        // Animation timer to keep text labels facing the camera.
        object : AnimationTimer() {
            override fun handle(now: Long) {
                fun updateText(node: Node) {
                    if (node is Text) {
                        node.transforms.setAll(
                            Rotate(-rotateY.angle, Rotate.Y_AXIS),
                            Rotate(-rotateX.angle, Rotate.X_AXIS)
                        )
                    }
                    if (node is Parent) node.childrenUnmodifiable.forEach { updateText(it) }
                }
                updateText(sceneGroup)
            }
        }.start()

        primaryStage.title = "Generalized Tetrahedron with Collapse Toggles"
        primaryStage.scene = scene
        primaryStage.show()
    }
}

fun main(args: Array<String>) {
    Application.launch(GraphTetrahedronApp::class.java, *args)
}
