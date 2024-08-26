package features.rtandim

class LinearEquation(private val m: Double, private val b: Double) {
    operator fun invoke(): (Double) -> Double {
        return { x -> m * x + b }
    }
}